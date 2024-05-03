rule baseRecalibrator:
    input:
        bam=wrkdir / "alignments" / "{sample}.cons.filtered.realigned.bam",
        dbsnp=dbsnp,
        genome=genome,
    output:
        table=wrkdir / "metrics" / "{sample}_recal_data.table",
        bam=temp(wrkdir / "alignments" / "{sample}_dedup.recall.bam"),
        analyse_covariates=wrkdir / "metrics" / "{sample}_covariates.pdf",
    conda:
        "../envs/gatk.yaml"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
    log:
        logdir / "gatk/{sample}_recal.log",
    message:
        "Recalibrating with GATK BaseRecalibrator"
    shell:
        " ( "
        " gatk BaseRecalibrator -I {input.bam} -R {input.genome} "
        " --known-sites {input.dbsnp} "
        " -O {output.table} && "
        " gatk ApplyBQSR -I {input.bam} -R {genome} --bqsr-recal-file {output.table} -O {output.bam} && "
        " gatk AnalyzeCovariates "
        " -bqsr {output.table} "
        " -plots {output.analyse_covariates} "
        " ) &> {log} "


rule sort_index:
    input:
        bam=wrkdir / "alignments" / "{sample}_dedup.recall.bam",
    output:
        bam=wrkdir / "alignments" / "{sample}_dedup.recall.sorted.bam",
        bai=wrkdir / "alignments" / "{sample}_dedup.recall.sorted.bam.bai",
    conda:
        "../envs/samtools.yaml"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=24 * 60,
        nodes=1,
    log:
        logdir / "samtools/{sample}_sort.log",
    message:
        "Sorting and indexing recalibrated bam file"
    shell:
        " ( "
        " samtools sort --threads 8 -o {output.bam}##idx##{output.bai} {input.bam} --write-index"
        " ) &> {log} "
