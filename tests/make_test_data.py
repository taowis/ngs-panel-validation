"""Build the small test data set in tests/data/ from the full HG001 inputs.

The pipeline calls every chromosome 1-22, X and Y, so the test reference keeps
those 24 contig names, each holding one window of hs37d5: a panel gene (or two
neighbouring genes) with 5 kb either side where the chromosome has one, else
30 kb starting at the capture interval nearest the middle of the chromosome. Coordinates in every BED
and VCF are shifted to the window. Reads are the pairs from the NIST7035 and
NIST7086 BAMs (pipeline v2) with both mates in a window, duplicates kept so
MarkDuplicates has work to do.

Needs the full data (scripts/00_download.sh panel, reference, reads) and the
results_v2 BAMs; the output is committed, so CI does not need any of them.

    python tests/make_test_data.py
"""
import gzip
import os
import subprocess

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
os.chdir(ROOT)
OUT = "tests/data"
REF = "data/reference/hs37d5.fa"
CAPTURE = "data/raw/nexterarapidcapture_expandedexome_targetedregions.bed.gz"
TRUTH_VCF = "data/truth/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED = "data/truth/HG001_GRCh37_1_22_v4.2.1_benchmark.bed"
TARGETS = "data/panels/panel_targets.bed"
BY_GENE = "data/panels/panel_targets_by_gene.bed"
EVALUABLE = "data/panels/panel_evaluable.bed"
BAMS = {"NIST7035": "results_v2/bam/NIST7035_run1.md.bam", "NIST7086": "results_v2/bam/NIST7086_run1.md.bam"}
CHROMS = [str(c) for c in range(1, 23)] + ["X", "Y"]
GENES = {"1": ["MUTYH"], "2": ["EPCAM", "MSH2"], "3": ["MLH1"], "5": ["APC"], "7": ["PMS2"], "10": ["PTEN"],
         "11": ["ATM"], "12": ["POLE"], "13": ["BRCA2"], "15": ["GREM1"], "16": ["PALB2"], "17": ["BRCA1"],
         "18": ["SMAD4"], "19": ["STK11"], "22": ["CHEK2"]}
PAD, DEFAULT = 5000, 30000


def sh(cmd):
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, capture_output=True, text=True).stdout


def bed(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(("track", "browser", "#")):
                continue
            f = line.rstrip("\n").split("\t")
            yield [f[0].removeprefix("chr"), int(f[1]), int(f[2])] + f[3:]


lengths = dict(l.split()[:2] for l in open(REF + ".fai"))
capture = list(bed(CAPTURE))
windows = {}
for c in CHROMS:
    if c in GENES:
        iv = [r for r in bed(BY_GENE) if r[0] == c and r[3] in GENES[c]]
        windows[c] = (min(r[1] for r in iv) - PAD, max(r[2] for r in iv) + PAD, ",".join(GENES[c]))
    else:
        mid = int(lengths[c]) // 2
        start = min((r[1] for r in capture if r[0] == c), key=lambda x: abs(x - mid))
        windows[c] = (start - 5000, start - 5000 + DEFAULT, "capture")

os.makedirs(f"{OUT}/ref", exist_ok=True)
with open(f"{OUT}/windows.tsv", "w") as fh:
    fh.write("contig\ths37d5_start\ths37d5_end\tcontent\n")
    for c, (s, e, what) in windows.items():
        fh.write(f"{c}\t{s}\t{e}\t{what}\n")
win_bed = f"{OUT}/windows.hs37d5.bed"
with open(win_bed, "w") as fh:
    for c, (s, e, _) in windows.items():
        fh.write(f"{c}\t{s}\t{e}\n")

# reference: one window per contig, renamed to the chromosome
with open(f"{OUT}/ref/hs37d5.fa", "w") as fh:
    for c, (s, e, _) in windows.items():
        seq = sh(f"samtools faidx {REF} {c}:{s + 1}-{e} | tail -n +2")
        fh.write(f">{c}\n{seq}")


def shift_bed(src, dst, gz=False):
    rows = []
    for r in bed(src):
        if r[0] not in windows:
            continue
        s, e, _ = windows[r[0]]
        a, b = max(r[1], s), min(r[2], e)
        if a < b:
            rows.append([r[0], a - s, b - s] + r[3:])
    rows.sort(key=lambda r: (CHROMS.index(r[0]), r[1]))
    text = "".join("\t".join(map(str, r)) + "\n" for r in rows)
    if gz:
        with gzip.open(dst, "wt") as fh:
            fh.write(text)
    else:
        open(dst, "w").write(text)


shift_bed(CAPTURE, f"{OUT}/capture.bed.gz", gz=True)
shift_bed(TRUTH_BED, f"{OUT}/HG001_test_truth.bed")
shift_bed(TARGETS, f"{OUT}/panel_targets.bed")
shift_bed(EVALUABLE, f"{OUT}/panel_evaluable.bed")

# truth VCF: records starting inside a window, positions shifted, contig lengths rewritten
body = sh(f"bcftools view --no-version -H -R {win_bed} {TRUTH_VCF}").splitlines()
header = [h for h in sh(f"bcftools view --no-version -h {TRUTH_VCF}").splitlines() if not h.startswith("##contig")]
with open(f"{OUT}/HG001_test_truth.vcf", "w") as fh:
    fh.write("\n".join(header[:-1]) + "\n")
    for c, (s, e, _) in windows.items():
        fh.write(f"##contig=<ID={c},length={e - s}>\n")
    fh.write(header[-1] + "\n")
    for line in body:
        f = line.split("\t")
        s, e, _ = windows[f[0]]
        pos = int(f[1])
        if s < pos and pos + len(f[3]) - 1 <= e:
            f[1] = str(pos - s)
            fh.write("\t".join(f) + "\n")
sh(f"bgzip -f {OUT}/HG001_test_truth.vcf && tabix -f -p vcf {OUT}/HG001_test_truth.vcf.gz")

# reads: pairs with both mates in a window
for lib, bam in BAMS.items():
    pre = f"{OUT}/{lib}"
    sh(f"samtools view -u -F 0x900 -M -L {win_bed} {bam} | samtools collate -u -O - "
       f"| samtools fastq -n -s /dev/null -0 /dev/null -1 {pre}_R1.fastq.gz -2 {pre}_R2.fastq.gz -")
os.remove(win_bed)
