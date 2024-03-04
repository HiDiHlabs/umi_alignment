rule fgbio:
    input:
        alignment=wrkdir / "alignments" / "{run_id}" / "{sample}_aln_{lane}.bam",
        fastq_r2=wrkdir / "fastq" / "{run_id}" / "{sample}_R2_{lane}.fastq.gz",
    output:
        temp(wrkdir / "alignments" / "{run_id}" / "{sample}_aln_{lane}_umi_annot.bam"),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=4 * 60,
        nodes=1,
    conda:
        "../envs/fgbio.yaml"
    log:
        logdir / "fgbio" / "{run_id}_{sample}_{lane}.log",
    message:
        "Annotating BAM with UMIs from fastq file."
    shell:
        "fgbio AnnotateBamWithUmis -i {input.alignment} -f {input.fastq_r2} -o {output} -t RX -q UQ -s true &> {log}"


rule group_reads:
    input:
        bam=wrkdir / "alignments" / "{sample}_mate_fix.bam",
    output:
        bam=temp(
            wrkdir / "alignments" / "{sample}_merged_aln_umi_annot_sorted_grouped.bam"
        ),
        stats=wrkdir / "metrics" / "{sample}.grouped-family-sizes.txt",
    params:
        allowed_edits=1,
    threads: 2
    resources:
        mem_mb=8000,
        runtime=24 * 60,
        nodes=1,
    conda:
        "../envs/fgbio.yaml"
    log:
        logdir / "fgbio/group_{sample}.log",
    message:
        "Grouping reads by UMI and position for consensus calling."
    shell:
        "fgbio -Xmx{resources.mem_mb}m --compression 1 --async-io GroupReadsByUmi "
        "--input {input.bam} "
        "--strategy Adjacency "
        "--edits {params.allowed_edits} "
        "--output {output.bam} "
        "--family-size-histogram {output.stats} "
        "&> {log} "


rule call_consensus_reads:
    input:
        bam=wrkdir / "alignments" / "{sample}_merged_aln_umi_annot_sorted_grouped.bam",
    output:
        bam=temp(wrkdir / "alignments" / "{sample}.cons.unmapped.bam"),
    params:
        min_reads=3,
        min_base_qual=20,
    threads: 4
    resources:
        mem_mb=4000,
        runtime=24 * 60,
        nodes=1,
    conda:
        "../envs/fgbio.yaml"
    log:
        logdir / "fgbio/call_consensus_reads.{sample}.log",
    message:
        "Calling consensus reads from grouped reads."
    shell:
        "fgbio -Xmx{resources.mem_mb}m --compression 0 CallMolecularConsensusReads "
        "--input {input.bam} "
        "--output {output.bam} "
        "--min-reads {params.min_reads} "
        "--min-input-base-quality {params.min_base_qual} "
        "--threads {threads} "
        "&> {log}"


rule filter_consensus_reads:
    input:
        bam=wrkdir / "alignments" / "{sample}.cons.unmapped.bam",
        ref=genome,
    output:
        bam=temp(wrkdir / "alignments" / "{sample}.cons.filtered.bam"),
    params:
        min_reads=3,
        min_base_qual=40,
        max_error_rate=0.2,
    threads: 8
    resources:
        mem_mb=8000,
        runtime=24 * 60,
        nodes=1,
    conda:
        "../envs/fgbio.yaml"
    log:
        logdir / "fgbio" / "filter_consensus_reads.{sample}.log",
    message:
        "Filtering consensus reads and sorting into coordinate order."
    shell:
        "("
        "samtools sort -n -u {input.bam} | "
        "fgbio -Xmx{resources.mem_mb}m --compression 1 FilterConsensusReads "
        "--input /dev/stdin "
        "--output {output.bam} "
        "--ref {input.ref} "
        "--min-reads {params.min_reads} "
        "--min-base-quality {params.min_base_qual} "
        "--max-base-error-rate {params.max_error_rate} "
        ") &> {log} "
