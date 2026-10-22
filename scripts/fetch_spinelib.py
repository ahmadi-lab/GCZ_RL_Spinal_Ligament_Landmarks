"""Download the SpineLib files of the released SLD implementation (needed only for SLD --fresh).

SpineLib belongs to the SpineToolkit 3D Slicer extension (https://github.com/VisSim-UniKO/SpineToolkit).
The repository has no licence file, so it is not redistributed here. This script fetches the
21 SpineLib files used by code/ligament/sld_official.py from the pinned commit, unmodified, and
verifies each against its git blob SHA-1 (scripts/spinelib_files.csv).

Usage:  python scripts/fetch_spinelib.py
"""
import csv
import hashlib
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
COMMIT = "62f39cffdc3957f2de35b5647c4e7bc9c79f9d28"
URL = f"https://raw.githubusercontent.com/VisSim-UniKO/SpineToolkit/{COMMIT}/"
DEST = REPO / "code" / "third_party"


def blob_sha1(data: bytes) -> str:
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def main() -> int:
    bad = 0
    rows = list(csv.DictReader(open(Path(__file__).with_name("spinelib_files.csv"), encoding="utf-8")))
    for r in rows:
        with urllib.request.urlopen(URL + r["path"], timeout=60) as resp:
            data = resp.read()
        ok = blob_sha1(data) == r["git_blob_sha1"]
        bad += not ok
        out = DEST / r["path"]
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(data)
        print(("ok  " if ok else "BAD ") + r["path"])
    print(f"{len(rows) - bad}/{len(rows)} SpineLib files verified (SpineToolkit@{COMMIT[:7]})")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
