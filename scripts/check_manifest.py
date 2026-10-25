"""Compare md5 and size of downloaded files with data/MANIFEST.tsv.

Files listed in the manifest but not on disk are skipped (not downloaded yet).
Exits non-zero on any mismatch. Run it after scripts/00_download.sh, and again
after re-downloading data that was deleted to free space, before rerunning an
analysis on it.

    python scripts/check_manifest.py
"""
import csv
import hashlib
import os
import sys

root = os.path.join(os.path.dirname(__file__), "..")
bad = 0
with open(os.path.join(root, "data", "MANIFEST.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        path = os.path.join(root, row["path"])
        if not os.path.exists(path):
            print(f"missing  {row['path']}")
            continue
        h = hashlib.md5()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        ok = h.hexdigest() == row["md5"] and os.path.getsize(path) == int(row["size_bytes"])
        bad += not ok
        print(f"{'ok' if ok else 'MISMATCH':8} {row['path']}")
sys.exit(1 if bad else 0)
