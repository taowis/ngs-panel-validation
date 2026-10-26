"""Compare a run on the test data with tests/expected/test_run.tsv.

Checked: no failed task in the trace, the VCF records of each sample (md5 of
the non-header lines, so paths and dates in the header do not matter), the
vcfeval counts for every sample, region, variant type and mode, and the
within-run comparison. Exits non-zero on any difference.

    python tests/check_test_run.py [results_test_data]
    python tests/check_test_run.py --write    # after an intended change, with the reason in the PR
"""
import argparse
import csv
import glob
import gzip
import hashlib
import os
import sys

p = argparse.ArgumentParser()
p.add_argument("outdir", nargs="?", default="results_test_data")
p.add_argument("--expected", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "expected/test_run.tsv"))
p.add_argument("--write", action="store_true")
args = p.parse_args()

trace = list(csv.DictReader(open(f"{args.outdir}/pipeline_info/trace.tsv"), delimiter="\t"))
failed = [r["name"] for r in trace if r["status"] not in ("COMPLETED", "CACHED")]
if failed:
    sys.exit(f"failed tasks: {', '.join(failed)}")

got = {}
for vcf in sorted(glob.glob(f"{args.outdir}/vcf/*.vcf.gz")):
    h = hashlib.md5()
    with gzip.open(vcf, "rb") as fh:
        for line in fh:
            if not line.startswith(b"#"):
                h.update(line)
    got[f"vcf_records_md5\t{os.path.basename(vcf)[:-7]}"] = h.hexdigest()
for f in sorted(glob.glob(f"{args.outdir}/vcfeval/*.summary.txt")):
    key = os.path.basename(f)[:-len(".summary.txt")]
    rows = [l.split() for l in open(f) if l.split() and l.split()[0] == "None"]
    got[f"vcfeval_tp_tp_fp_fn\t{key}"] = ",".join(rows[0][1:5]) if rows else "no variants"
for f in sorted(glob.glob(f"{args.outdir}/reproducibility/within_run.*.tsv")):
    r = next(csv.DictReader(open(f), delimiter="\t"))
    got[f"within_run_differing_lines\t{r['run_a']}_vs_{r['run_b']}"] = r["differing_lines"]

if args.write:
    os.makedirs(os.path.dirname(args.expected), exist_ok=True)
    with open(args.expected, "w") as fh:
        fh.write("check\titem\tvalue\n")
        for k, v in got.items():
            fh.write(f"{k}\t{v}\n")
    print(f"wrote {len(got)} values to {args.expected}")
    sys.exit()

expected = {}
for line in list(open(args.expected))[1:]:
    check, item, value = line.rstrip("\n").split("\t")
    expected[f"{check}\t{item}"] = value
bad = 0
for k in sorted(set(expected) | set(got)):
    if expected.get(k) != got.get(k):
        bad += 1
        print(f"DIFFERS  {k.replace(chr(9), ' ')}: expected {expected.get(k)}, got {got.get(k)}")
print(f"{len(expected) - bad} of {len(expected)} checks match" if not bad else f"{bad} differences")
sys.exit(1 if bad else 0)
