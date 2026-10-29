"""Compare the outputs of a run with the committed versions.

Run it in a clone of the repository after RUN_ALL.ps1 or run_all.sh:

    python scripts/compare_outputs.py

It lists every committed file that the run changed and classifies the change:

  same data        CSV and JSON values equal within 1e-9 (last-digit floating-point noise)
  rendering noise  PNG/JPG images and PDF pages (rasterized with pdftoppm at 100 dpi) in which
                   fewer than 0.2 % of the pixels differ by more than 24/255
  provenance       run logs, run times and dates, SVG files (their PDF and PNG twins are compared)
  CHANGED          anything else

Exit code 0 if nothing is CHANGED, 1 otherwise.
"""
import io
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
TOL = 1e-9
# JSON keys that record when or where a run happened, not what it computed
PROVENANCE_KEYS = re.compile(r"time|date|elapsed|duration|seconds|started|finished|host|platform|"
                             r"python|version|gcz_geometry|status|log", re.I)
PDFTOPPM = shutil.which("pdftoppm")


def git(*args, text=True):
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=text).stdout


def committed(f):
    return git("show", f"HEAD:{f}", text=False)


def pixels(a, b):
    """Share of pixels that differ by more than 24/255 in any channel, or None if the sizes differ."""
    if a.shape != b.shape:
        return None
    return float((np.abs(a - b).max(axis=2) > 24).mean())


def image(data):
    return np.asarray(Image.open(io.BytesIO(data)).convert("RGB"), float)


def pages(data):
    with tempfile.TemporaryDirectory() as t:
        (Path(t) / "a.pdf").write_bytes(data)
        subprocess.run([PDFTOPPM, "-r", "100", "-png", str(Path(t) / "a.pdf"), str(Path(t) / "p")], check=True)
        return [np.asarray(Image.open(p).convert("RGB"), float) for p in sorted(Path(t).glob("p*.png"))]


def json_diff(a, b, path=""):
    """(largest numeric difference, list of differing non-provenance paths, list of provenance paths)."""
    if isinstance(a, dict) and isinstance(b, dict):
        big, bad, prov = 0.0, [], []
        for k in set(a) | set(b):
            p = f"{path}/{k}"
            if k not in a or k not in b:
                (prov if PROVENANCE_KEYS.search(str(k)) else bad).append(p)
                continue
            if PROVENANCE_KEYS.search(str(k)) and a[k] != b[k]:
                prov.append(p)
                continue
            m, x, y = json_diff(a[k], b[k], p)
            big, bad, prov = max(big, m), bad + x, prov + y
        return big, bad, prov
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return 0.0, [f"{path} (length {len(a)} -> {len(b)})"], []
        big, bad, prov = 0.0, [], []
        for i, (x, y) in enumerate(zip(a, b)):
            m, u, v = json_diff(x, y, f"{path}[{i}]")
            big, bad, prov = max(big, m), bad + u, prov + v
        return big, bad, prov
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if math.isnan(a) and math.isnan(b):
            return 0.0, [], []
        d = abs(a - b)
        return d, ([path] if d > TOL else []), []
    return 0.0, ([] if a == b else [path]), []


def csv_diff(a, b):
    """Largest numeric difference, or a reason why the tables differ."""
    x, y = pd.read_csv(io.BytesIO(a)), pd.read_csv(io.BytesIO(b))
    if list(x.columns) != list(y.columns) or x.shape != y.shape:
        return None, f"shape or columns changed {x.shape} -> {y.shape}"
    big = 0.0
    for c in x.columns:
        u, v = x[c], y[c]
        if pd.api.types.is_numeric_dtype(u) and pd.api.types.is_numeric_dtype(v):
            u, v = u.to_numpy(float), v.to_numpy(float)
            if not np.array_equal(np.isnan(u), np.isnan(v)):
                return None, f"column {c}: missing values changed"
            d = np.abs(u - v)[~np.isnan(u)]
            big = max(big, float(d.max()) if d.size else 0.0)
        elif not u.astype(str).equals(v.astype(str)):
            return None, f"column {c}: text changed"
    return big, None


def classify(f):
    p = ROOT / f
    old, new = committed(f), p.read_bytes()
    suf = p.suffix.lower()
    if suf == ".log" or "/logs/" in f"/{f}":
        return "provenance", "run log"
    if "log" in p.name.lower() and suf == ".txt":
        stamp = re.compile(rb"^\[\d\d:\d\d:\d\d\] ?", re.M)
        if stamp.sub(b"", old).replace(b"\r\n", b"\n") == stamp.sub(b"", new).replace(b"\r\n", b"\n"):
            return "provenance", "run log, time stamps only"
    if suf == ".svg":
        return "provenance", "SVG with new embedded dates (compared through its PDF and PNG twins)"
    if suf in (".png", ".jpg", ".jpeg"):
        fr = pixels(image(old), image(new))
        if fr is None:
            return "CHANGED", "image size changed"
        return ("rendering noise" if fr < 0.002 else "CHANGED"), f"{100 * fr:.3f} % of pixels differ"
    if suf == ".pdf":
        if not PDFTOPPM:
            return "provenance", "PDF not rasterized (pdftoppm not found); compare its PNG twin or the page images"
        a, b = pages(old), pages(new)
        if len(a) != len(b):
            return "CHANGED", f"{len(a)} -> {len(b)} pages"
        fr = [pixels(x, y) for x, y in zip(a, b)]
        if any(x is None for x in fr):
            return "CHANGED", "page size changed"
        return ("rendering noise" if max(fr) < 0.002 else "CHANGED"), f"{len(a)} pages, {100 * max(fr):.3f} % of pixels differ"
    if suf == ".csv":
        big, why = csv_diff(old, new)
        if why:
            return "CHANGED", why
        return ("same data" if big <= TOL else "CHANGED"), f"largest difference {big:.1e}"
    if suf == ".json":
        big, bad, prov = json_diff(json.loads(old), json.loads(new))
        if bad:
            return "CHANGED", "differs at " + ", ".join(bad[:5]) + (" ..." if len(bad) > 5 else "")
        if prov and big <= TOL:
            return "provenance", "only " + ", ".join(sorted(set(re.sub(r"\[\d+\]", "[]", q) for q in prov))[:4])
        return "same data", f"largest difference {big:.1e}"
    if old.replace(b"\r\n", b"\n") == new.replace(b"\r\n", b"\n"):
        return "same data", "line endings only"
    return "CHANGED", (git("diff", "--stat", "--", f).strip().splitlines() or ["text changed"])[-1].strip()


def main():
    status = git("status", "--porcelain", "-z").split("\0")
    changed = [s[3:] for s in status if s[:2] in (" M", "M ", "MM")]
    deleted = [s[3:] for s in status if s[:2] in (" D", "D ")]
    groups = {"CHANGED": [], "same data": [], "rendering noise": [], "provenance": []}
    for f in changed:
        kind, note = classify(f)
        groups[kind].append(f"{f}: {note}")
    groups["CHANGED"] += [f"{f}: deleted" for f in deleted]
    print(f"committed files changed by the run: {len(changed) + len(deleted)}")
    for kind in ("CHANGED", "same data", "rendering noise", "provenance"):
        print(f"\n{kind}: {len(groups[kind])}")
        for line in groups[kind]:
            print("  " + line)
    ok = not groups["CHANGED"]
    print("\nRESULT: " + ("every output equals the committed version (within the tolerances above)" if ok
                          else "some outputs differ from the committed versions (see CHANGED)"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
