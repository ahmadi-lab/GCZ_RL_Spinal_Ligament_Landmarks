"""Execute the original GCZ notebook on the packaged data and compare with its saved outputs.

The notebook `GCZ_journal_FINAL_57_20_37_eval9_ORIGINAL.ipynb` is the source of dissertation
Table 4.2 and of journal-paper Tables 1-4 (dataset, alignment, matching, per-ligament
matching). It is left unchanged. This script runs its code cells in order with only these
run-time substitutions:

  * BASE / OUTPUT_BASE point to ../../dataset/raw/VerSe_Ligaments_GT and ./outputs
    (the original absolute paths no longer exist);
  * figures are written as HTML files but no browser window is opened;
  * the IPython line magic `%matplotlib inline` is skipped and `display()` is a no-op.

Outputs (./outputs):
  GCZ_notebook_executed.ipynb   the notebook with the new outputs
  cell_outputs/cell_NN.txt      text output of every code cell
  comparison_with_saved_outputs.csv   per cell: identical / differs vs. the saved outputs
  tables/                       the published tables extracted from the regenerated outputs

Usage:  python run_original_notebook.py      (about 30-60 min)
"""
from __future__ import annotations

import contextlib
import copy
import io
import json
import os
import re
import shutil
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NB = HERE / "GCZ_journal_FINAL_57_20_37_eval9_ORIGINAL.ipynb"
DATA = Path(os.environ.get("GCZRL_RAW_ROOT", ROOT / "dataset" / "raw")) / "VerSe_Ligaments_GT"
OUT = HERE / "outputs"
CELL_OUT = OUT / "cell_outputs"


def patched_source(i: int, src: str) -> str:
    if i == 1:
        src = re.sub(r'^BASE = r".*"$', lambda m: f'BASE = r"{DATA}"', src, count=1, flags=re.M)
        src = re.sub(r'^OUTPUT_BASE = r".*"$', lambda m: f'OUTPUT_BASE = r"{OUT}"', src, count=1, flags=re.M)
        src = src.replace('webbrowser.open("file://" + os.path.abspath(fname))', "pass")
    return "\n".join("" if ln.lstrip().startswith(("%", "!")) else ln for ln in src.splitlines())


def saved_text(cell) -> str:
    return "".join("".join(o.get("text", [])) for o in cell.get("outputs", []) if o.get("output_type") == "stream")


def normalise(s: str) -> str:
    """Remove run-specific paths and timings before comparing."""
    s = s.replace(str(DATA), "<BASE>").replace(str(OUT), "<OUT>")
    s = re.sub(r"C:\\Users\\mohse\\Documents\\ligament_AI - Copy \(12\)\\mohsen\\Dataset\\VerSe_Ligaments_GT", "<BASE>", s)
    s = re.sub(r"C:\\Users\\mohse\\Documents\\ligament_AI - Copy \(12\)\\mohsen", "<OUT>", s)
    s = re.sub(r"\d+(\.\d+)?\s*s\b", "<t>s", s)
    return s.strip()


def main():
    import matplotlib
    matplotlib.use("Agg")
    import pandas as pd

    OUT.mkdir(exist_ok=True)
    CELL_OUT.mkdir(parents=True, exist_ok=True)
    # the split used by the dissertation, read explicitly instead of being rebuilt
    shutil.copy(ROOT / "dataset" / "split_manifest_57.csv", OUT / "gcz_57_split_manifest.csv")
    before = {p for p in DATA.parent.rglob("*") if p.is_file()}

    nb = json.loads(NB.read_text(encoding="utf-8"))
    run_nb = copy.deepcopy(nb)
    ns = {"__name__": "__main__", "display": lambda *a, **k: None}
    os.chdir(OUT)
    rows = []
    for i, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != "code":
            continue
        buf, t0, err = io.StringIO(), time.time(), ""
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            try:
                exec(compile(patched_source(i, "".join(cell["source"])), f"<cell {i}>", "exec"), ns)
            except Exception:
                err = traceback.format_exc()
        text = buf.getvalue() + (("\n" + err) if err else "")
        import matplotlib.pyplot as plt
        plt.close("all")
        (CELL_OUT / f"cell_{i:02d}.txt").write_text(text, encoding="utf-8")
        run_nb["cells"][i]["outputs"] = [{"output_type": "stream", "name": "stdout", "text": text.splitlines(True)}]
        old = saved_text(cell)
        same = normalise(old) == normalise(text)
        rows.append({"cell": i, "first_line": "".join(cell["source"]).strip().split("\n")[0][:70],
                     "seconds": round(time.time() - t0, 1), "error": bool(err),
                     "saved_output_chars": len(old), "new_output_chars": len(text),
                     "identical_to_saved_output": same})
        print(f"cell {i:2d}  {rows[-1]['seconds']:7.1f} s  {'ERROR' if err else ('identical' if same else 'differs')}",
              file=sys.__stdout__, flush=True)

    (OUT / "GCZ_notebook_executed.ipynb").write_text(json.dumps(run_nb, indent=1), encoding="utf-8")
    pd.DataFrame(rows).to_csv(OUT / "comparison_with_saved_outputs.csv", index=False)

    # files the notebook writes next to the meshes are moved out of the dataset folder
    moved = OUT / "files_written_by_notebook_into_dataset"
    for p in sorted({p for p in DATA.parent.rglob("*") if p.is_file()} - before):
        dst = moved / p.relative_to(DATA.parent)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(p), dst)
    for d in sorted((p for p in DATA.rglob("*") if p.is_dir()), reverse=True):
        if not any(d.iterdir()):
            d.rmdir()
    n_same = sum(r["identical_to_saved_output"] for r in rows)
    print(f"\n{n_same}/{len(rows)} code cells give output identical to the saved notebook; "
          f"errors: {sum(r['error'] for r in rows)}", file=sys.__stdout__)


if __name__ == "__main__":
    main()
