nextflow.enable.dsl = 2

// Germline exome pipeline (BWA-MEM, MarkDuplicates, HaplotypeCaller, hard
// filters) plus the validation steps: vcfeval against GIAB, coverage over the
// panel, within-run and between-library comparisons.

def CHROMS = (1..22).collect { it.toString() } + ['X', 'Y']

process PREP_REGIONS {
    input:
    path capture
    path truth_bed

    output:
    path 'capture.bed', emit: capture
    path 'exome_evaluable.bed', emit: exome

    script:
    """
    gzip -dc ${capture} | sed 's/^chr//' | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > capture.bed
    bedtools intersect -a capture.bed -b ${truth_bed} | sort -k1,1 -k2,2n | bedtools merge > exome_evaluable.bed
    """
}

process BWA_MEM {
    tag "${sample}.${lane}"
    cpus 6
    memory '7 GB'

    input:
    tuple val(sample), val(library), val(lane), path(r1), path(r2)
    path ref

    output:
    tuple val(sample), path("${sample}.${lane}.bam")

    script:
    // -K fixes the batch size so output does not depend on thread count.
    // bwa holds about 5.5 GB of index; sort uses 2 x 768 MB on top of that.
    """
    bwa mem -t ${task.cpus} -K 100000000 -Y \\
        -R '@RG\\tID:${library}.${lane}\\tSM:${params.genome}\\tLB:${library}\\tPL:ILLUMINA\\tPU:${library}.${lane}' \\
        hs37d5.fa ${r1} ${r2} \\
      | samtools sort -@ 1 -m 768M -o ${sample}.${lane}.bam -
    """
}

process MARKDUP {
    tag "${sample}"
    cpus 2
    memory '5 GB'
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.txt'
    // hard links, so the BAMs take no extra space next to the work directory
    publishDir "${params.outdir}/bam", mode: 'link', pattern: '*.md.bam*'

    input:
    tuple val(sample), path(bams)

    output:
    tuple val(sample), path("${sample}.md.bam"), path("${sample}.md.bam.bai"), emit: bam
    path "${sample}.markdup_metrics.txt", emit: metrics

    script:
    def inputs = bams.collect { "-I ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx4g" MarkDuplicates ${inputs} \\
        -O ${sample}.md.bam -M ${sample}.markdup_metrics.txt
    samtools index ${sample}.md.bam
    """
}

process MEAN_DEPTH {
    tag "${sample}"
    cpus 4
    memory '4 GB'

    input:
    tuple val(sample), path(bam), path(bai)
    path capture

    output:
    tuple val(sample), path(bam), path(bai), env(MEAN)

    script:
    // same filters as COVERAGE: MAPQ >= 20, duplicates excluded
    """
    mosdepth -t ${task.cpus} -n -Q 20 --by ${capture} full ${bam}
    MEAN=\$(awk '\$1 == "total_region" { print \$4 }' full.mosdepth.summary.txt)
    """
}

process DOWNSAMPLE {
    tag "${sample}.${depth}x"
    cpus 2
    memory '2 GB'

    input:
    tuple val(sample), path(bam), path(bai), val(mean), val(depth)

    output:
    tuple val("${sample}_${depth}x"), path("${sample}_${depth}x.bam"), path("${sample}_${depth}x.bam.bai")

    script:
    // samtools -s takes seed.fraction; read pairs are kept or dropped together.
    // The seed is fixed so a rerun keeps the same reads; the fraction is the
    // target depth over the mean measured by MEAN_DEPTH.
    def s = '42' + String.format('%.6f', depth / mean.toDouble()).substring(1)
    """
    samtools view -@ ${task.cpus} -b -s ${s} ${bam} -o ${sample}_${depth}x.bam
    samtools index ${sample}_${depth}x.bam
    """
}

process HAPLOTYPECALLER {
    tag "${sample}.${chrom}"
    cpus 1
    memory '3 GB'
    maxForks 3

    input:
    tuple val(sample), path(bam), path(bai), val(chrom)
    path ref
    path capture

    output:
    tuple val(sample), val(chrom), path("${sample}.${chrom}.vcf.gz")

    script:
    // one chromosome per task; maxForks 3 x 3 GB keeps calling inside the 9 GB
    // executor limit. -ip 100 lets reads at the edge of a capture interval be
    // assembled with their neighbours.
    """
    awk -v c=${chrom} '\$1 == c' ${capture} > regions.bed
    gatk --java-options "-Xmx2g" HaplotypeCaller -R hs37d5.fa -I ${bam} \\
        -L regions.bed -ip 100 -O ${sample}.${chrom}.vcf.gz
    """
}

process GATHER_FILTER {
    tag "${sample}"
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    tuple val(sample), val(chroms), path(vcfs)

    output:
    tuple val(sample), path("${sample}.vcf.gz"), path("${sample}.vcf.gz.tbi")

    script:
    // concatenate in reference order, then GATK's generic hard filters, except
    // that SNPs are not filtered on SOR (change 02, docs/reports/02-snv-filter-change).
    // Thresholds are the values in GATK's "Hard-filtering germline short variants"
    // article. Hard filters rather than VQSR, which needs more variants than a
    // single exome gives.
    def order = CHROMS.collect { c -> "${sample}.${c}.vcf.gz" }.join(' ')
    """
    bcftools concat --no-version ${order} -Oz -o raw.vcf.gz
    bcftools filter --no-version -m+ -s SNP_HARD \\
        -e 'TYPE="snp" && (INFO/QD<2 || INFO/FS>60 || INFO/MQ<40 || INFO/MQRankSum<-12.5 || INFO/ReadPosRankSum<-8)' \\
        raw.vcf.gz \\
      | bcftools filter --no-version -m+ -s INDEL_HARD \\
        -e 'TYPE!="snp" && (INFO/QD<2 || INFO/FS>200 || INFO/ReadPosRankSum<-20 || INFO/SOR>10)' \\
        -Oz -o ${sample}.vcf.gz
    tabix -p vcf ${sample}.vcf.gz
    """
}

process VCFEVAL {
    tag "${sample}.${region}.${vtype}.${mode}"
    cpus 2
    memory '4 GB'
    publishDir "${params.outdir}/vcfeval", mode: 'copy', pattern: '*.summary.txt'

    input:
    tuple val(sample), path(vcf), path(tbi), val(region), path(region_bed), val(vtype), val(mode)
    path truth_vcf
    path truth_tbi
    path truth_bed
    path sdf

    output:
    path "${sample}.${region}.${vtype}.${mode}.summary.txt"

    script:
    // genotype mode needs matching genotypes; allele mode (--squash-ploidy)
    // only matching alleles, so genotype concordance = TP_genotype / TP_allele
    def sel = vtype == 'snv' ? 'snps' : 'indels'
    def squash = mode == 'allele' ? '--squash-ploidy' : ''
    """
    bcftools view --no-version -v ${sel} ${truth_vcf} -Oz -o truth.vcf.gz && tabix truth.vcf.gz
    bcftools view --no-version -v ${sel} ${vcf} -Oz -o calls.vcf.gz && tabix calls.vcf.gz
    rtg vcfeval -b truth.vcf.gz -c calls.vcf.gz -t ${sdf} \\
        --bed-regions ${region_bed} --evaluation-regions ${truth_bed} \\
        ${squash} --threads ${task.cpus} -o out
    cp out/summary.txt ${sample}.${region}.${vtype}.${mode}.summary.txt
    """
}

process COVERAGE {
    tag "${sample}"
    cpus 4
    memory '4 GB'
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path targets
    path capture

    output:
    tuple val(sample), path("${sample}.panel_depth.tsv.gz"), emit: panel
    tuple val(sample), path("${sample}.callable.bed"), emit: callable
    path "${sample}.mosdepth.summary.txt", emit: summary

    script:
    // panel: per base, MAPQ >= 20 and base quality >= 20, duplicates excluded
    // exome: mosdepth regions at >= min_depth (MAPQ >= 20) for between-library comparison
    """
    samtools depth -a -b ${targets} -q 20 -Q 20 ${bam} | gzip > ${sample}.panel_depth.tsv.gz
    MOSDEPTH_Q0=LOW MOSDEPTH_Q1=CALLABLE mosdepth -t ${task.cpus} -n -Q 20 \\
        --by ${capture} --quantize 0:${params.min_depth}: ${sample} ${bam}
    gzip -dc ${sample}.quantized.bed.gz | awk '\$4 == "CALLABLE"' | cut -f1-3 > ${sample}.callable.bed
    """
}

process WITHIN_RUN {
    publishDir "${params.outdir}/reproducibility", mode: 'copy'

    input:
    tuple val(a), path(vcf_a), path(tbi_a), val(b), path(vcf_b), path(tbi_b)

    output:
    path "within_run.${a}_vs_${b}.tsv"

    script:
    """
    bcftools view -H ${vcf_a} > a.txt
    bcftools view -H ${vcf_b} > b.txt
    n_a=\$(wc -l < a.txt); n_b=\$(wc -l < b.txt)
    n_diff=\$(diff a.txt b.txt | grep -c '^[<>]' || true)
    printf 'run_a\\trun_b\\trecords_a\\trecords_b\\tdiffering_lines\\tmd5_a\\tmd5_b\\n' > within_run.${a}_vs_${b}.tsv
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' ${a} ${b} \$n_a \$n_b \$n_diff \$(md5sum a.txt | cut -c1-32) \$(md5sum b.txt | cut -c1-32) >> within_run.${a}_vs_${b}.tsv
    """
}

process BETWEEN_LIBRARY {
    cpus 2
    memory '4 GB'
    publishDir "${params.outdir}/reproducibility", mode: 'copy'

    input:
    tuple val(a), path(vcf_a), path(tbi_a), path(callable_a), val(b), path(vcf_b), path(tbi_b), path(callable_b)
    path targets
    path sdf

    output:
    path "between_library.${a}_vs_${b}.*.summary.txt"

    script:
    // both libraries at >= min_depth; exome-wide and inside the panel targets
    """
    bedtools intersect -a ${callable_a} -b ${callable_b} | sort -k1,1 -k2,2n | bedtools merge > both.bed
    bedtools intersect -a both.bed -b ${targets} | sort -k1,1 -k2,2n | bedtools merge > both_panel.bed
    for r in exome panel; do
        bed=both.bed; [ \$r = panel ] && bed=both_panel.bed
        out=between_library.${a}_vs_${b}.\$r.summary.txt
        if [ -s \$bed ]; then
            rtg vcfeval -b ${vcf_a} -c ${vcf_b} -t ${sdf} --bed-regions \$bed --threads ${task.cpus} -o out_\$r
            cp out_\$r/summary.txt \$out
        else
            echo "no bases at >= ${params.min_depth}x in both libraries" > \$out
        fi
    done
    """
}

workflow {
    ref = Channel.value(files(params.ref + '{,.fai,.amb,.ann,.bwt,.pac,.sa}') + [file(params.ref.replaceAll(/\.fa$/, '.dict'))])
    sdf = Channel.value(file(params.sdf, type: 'dir'))
    truth_vcf = Channel.value(file(params.truth_vcf))
    truth_tbi = Channel.value(file(params.truth_vcf + '.tbi'))
    truth_bed = Channel.value(file(params.truth_bed))
    targets = Channel.value(file(params.targets))
    panel_eval = Channel.value(file(params.evaluable))

    regions = PREP_REGIONS(file(params.capture), file(params.truth_bed))

    reads = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { r -> tuple(r.sample, r.library, r.lane, file(r.fastq_1), file(r.fastq_2)) }

    lanes = BWA_MEM(reads, ref).groupTuple()
    bams = MARKDUP(lanes).bam
    if (params.downsample) {
        // also call at lower mean depths over the capture (report 03)
        levels = Channel.fromList(params.downsample.toString().tokenize(',')*.toInteger())
        full = MEAN_DEPTH(bams, regions.capture)
        bams = bams.mix(DOWNSAMPLE(full.combine(levels).filter { it[4] < it[3].toDouble() }))
    }

    calls = HAPLOTYPECALLER(bams.combine(Channel.fromList(CHROMS)), ref, regions.capture)
        .map { s, c, v -> tuple(s, c, v) }
        .groupTuple(size: CHROMS.size())
    vcfs = GATHER_FILTER(calls)

    region_beds = Channel.of('panel').combine(panel_eval)
        .mix(Channel.of('exome').combine(regions.exome))
    jobs = vcfs.combine(region_beds).combine(Channel.of('snv', 'indel')).combine(Channel.of('genotype', 'allele'))
    VCFEVAL(jobs, truth_vcf, truth_tbi, truth_bed, sdf)

    cov = COVERAGE(bams, targets, regions.capture)

    by_name = vcfs.map { s, v, t -> tuple(s, v, t) }
    run1 = by_name.filter { it[0] == params.rep_a }
    run2 = by_name.filter { it[0] == params.rep_b }
    WITHIN_RUN(run1.combine(run2))

    with_cov = by_name.join(cov.callable)
    lib_a = with_cov.filter { it[0] == params.rep_a }
    lib_b = with_cov.filter { it[0] == params.lib_b }
    BETWEEN_LIBRARY(lib_a.combine(lib_b), targets, sdf)
}
