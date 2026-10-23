"""Report 03, analysis 7: sample identity, sex and trio consistency.

identity  each sample's PASS SNVs against every genome's truth set
          (vcfeval genotype mode, chromosomes 20 to 22 to keep it quick), over
          the sample's capture & that truth set's high-confidence regions
sex       heterozygous fraction of PASS SNVs on chrX outside the PARs, and
          chrY mean depth relative to the whole capture (mosdepth)
trio      HG002 (child), HG003, HG004: PASS genotypes at positions at >= 20x
          in all three; a site not called in a sample counts as homozygous
          reference there

    python scripts/09_identity.py
"""
import argparse
import os

from genomes import GENOMES, ROOT, load, sh, write_tsv

p = argparse.ArgumentParser()
p.add_argument("--report", default="docs/reports/03-more-reference-samples")
args = p.parse_args()
os.chdir(ROOT)
SEX = {"HG001": "female", "HG002": "male", "HG003": "male", "HG004": "female", "HG005": "male"}
PAR = "X:60001-2699520,X:154931044-155260560"  # GRCh37 pseudo-autosomal regions

gs = {n: load(n) for n in GENOMES}
samples = [(n, s) for n in GENOMES for s in gs[n]["samples"] if not s.endswith("_run2")]

# identity
id_rows = []
for gname, s in samples:
    g = gs[gname]
    work = f"{g['outdir']}/identity"
    os.makedirs(work, exist_ok=True)
    vcf = f"{g['outdir']}/vcf/{s}.vcf.gz"
    calls = f"{work}/{s}.snv.vcf.gz"
    if not os.path.exists(calls + ".tbi"):
        sh(f"bcftools view --no-version -f PASS -v snps -r 20,21,22 {vcf} -Oz -o {calls} && tabix -f {calls}")
    for tname in GENOMES:
        t = gs[tname]
        region = f"{work}/{s}.vs.{tname}.bed"
        out = f"{work}/{s}.vs.{tname}"
        if not os.path.exists(f"{out}/summary.txt"):
            sh(f"gzip -dc {g['capture']} | sed 's/^chr//' | cut -f1-3 | awk '$1 == 20 || $1 == 21 || $1 == 22' "
               f"| sort -k1,1 -k2,2n | bedtools merge | bedtools intersect -a - -b {t['truth_bed']} > {region}")
            sh(f"rm -rf {out} && bcftools view --no-version -v snps -r 20,21,22 {t['truth_vcf']} "
               f"-Oz -o {out}.truth.vcf.gz && tabix -f {out}.truth.vcf.gz "
               f"&& rtg vcfeval -b {out}.truth.vcf.gz -c {calls} -t data/reference/hs37d5.sdf "
               f"--bed-regions {region} --threads 4 -o {out}")
        for line in open(f"{out}/summary.txt"):
            f = line.split()
            if f and f[0] == "None":
                tpb, tpc, fp, fn = map(int, f[1:5])
                fm = f[7]
        id_rows.append([gname, s, tname, tpb, fp, fn, fm])
write_tsv(f"{args.report}/identity.tsv", ["genome", "sample", "truth_genome", "tp_baseline", "fp", "fn",
                                          "f_measure"], id_rows)

# sex
sex_rows = []
for gname, s in samples:
    g = gs[gname]
    gts = sh(f"bcftools view -f PASS -v snps -r X -t ^{PAR} {g['outdir']}/vcf/{s}.vcf.gz "
             f"| bcftools query -f '[%GT]\\n'").split()
    het = sum(1 for x in gts if x in ("0/1", "0|1", "1|0"))
    mean = {}
    for line in open(f"{g['outdir']}/coverage/{s}.mosdepth.summary.txt"):
        f = line.split()
        mean[f[0]] = f[3]
    y = float(mean.get("Y_region", 0)) / float(mean["total_region"])
    inferred = "male" if het / len(gts) < 0.2 and y > 0.1 else "female" if het / len(gts) > 0.3 and y < 0.05 \
        else "unclear"
    sex_rows.append([gname, s, len(gts), het, f"{het / len(gts):.3f}", f"{y:.3f}", inferred, SEX[gname],
                     "match" if inferred == SEX[gname] else "MISMATCH"])
write_tsv(f"{args.report}/sex_check.tsv", ["genome", "sample", "chrx_snvs", "chrx_het", "chrx_het_fraction",
                                           "chry_relative_depth", "inferred", "expected", "result"], sex_rows)

# trio
trio = ["HG002", "HG003", "HG004"]
work = f"{gs['HG002']['outdir']}/trio"
os.makedirs(work, exist_ok=True)
callable_bed = f"{work}/callable_all3.bed"
beds = [f"{gs[n]['outdir']}/coverage/{n}.callable.bed" for n in trio]
sh(f"bedtools intersect -a {beds[0]} -b {beds[1]} | bedtools intersect -a - -b {beds[2]} "
   f"| awk '$1 ~ /^[0-9]+$/' | sort -k1,1 -k2,2n | bedtools merge > {callable_bed}")
parts = []
for n in trio:
    out = f"{work}/{n}.pass.vcf.gz"
    sh(f"bcftools view --no-version -f PASS -R {callable_bed} {gs[n]['outdir']}/vcf/{n}.vcf.gz -Oz -o {out} "
       f"&& tabix -f {out}")
    parts.append(out)
merged = sh(f"bcftools merge --no-version -m all -0 {' '.join(parts)} "
            f"| bcftools query -f '%TYPE\\t[%GT\\t]\\n'").splitlines()
counts = {}
for line in merged:
    typ, *g = line.rstrip("\t").split("\t")
    vtype = "snv" if typ == "SNP" else "indel"
    al = [set(x.replace("|", "/").split("/")) for x in g]
    if any("." in a for a in al) or all(a == {"0"} for a in al):
        continue
    child = g[0].replace("|", "/").split("/")
    ok = (child[0] in al[1] and child[1] in al[2]) or (child[1] in al[1] and child[0] in al[2])
    n, e = counts.get(vtype, (0, 0))
    counts[vtype] = (n + 1, e + (not ok))
bp = sum(int(e) - int(b) for _, b, e in (l.split() for l in open(callable_bed)))
trio_rows = [["HG002", "HG003", "HG004", t, bp, n, e, f"{e / n:.4f}"] for t, (n, e) in sorted(counts.items())]
write_tsv(f"{args.report}/trio_mendelian.tsv", ["child", "father", "mother", "type", "callable_bp", "variant_sites",
                                               "inconsistent", "rate"], trio_rows)
for r in id_rows + sex_rows + trio_rows:
    print("\t".join(map(str, r)))
