


rule fix_mate:
    """
    Fixing mate information if required
    For some reason for some reads the mate information is not properly set.
    This can cause problems in downstream analysis.
    """
    input:
        bam=(
            wrkdir / "alignments" / "{sample}_merged_umi_annot.bam"
            if not read_structure
            else wrkdir / "alignments" / "{sample}_merged_umi_annot.bam"
        ),
    output:
        bam=temp(wrkdir / "alignments" / "{sample}_mate_fix.bam"),
    threads: 2
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/fgbio.yaml"
    log:
        logdir / "fgbio/fixmate_{sample}.log",
    message:
        "Fixing mate information if required"
    shell:
        "fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m --compression 1 --async-io SetMateInformation "
        "--input {input.bam} "
        "--output {output.bam} "
        "--allow-missing-mates true"
        "&> {log} "


rule group_reads:
    input:
        bam=wrkdir / "alignments" / "{sample}_mate_fix.bam",
    output:
        bam=temp(
            wrkdir / "alignments" / "{sample}_merged_aln_umi_annot_sorted_grouped.bam"
        ),
        stats=wrkdir / "metrics" / "{sample}.grouped-family-sizes.txt",
    params:
        strategy=group_strategy,
        allowed_edits=group_allowed_edits,
        min_mapq=group_min_mapq,
    threads: 2
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/fgbio.yaml"
    log:
        logdir / "fgbio/group_{sample}.log",
    message:
        "Grouping reads by UMI and position for consensus calling."
    shell:
        "fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m --compression 1 --async-io GroupReadsByUmi "
        "--input {input.bam} "
        "--strategy {params.strategy} "
        "--edits {params.allowed_edits} "
        "--min-map-q {params.min_mapq} "
        "--output {output.bam} "
        "--family-size-histogram {output.stats} "
        "&> {log} "


rule call_consensus_reads:
    input:
        bam=wrkdir / "alignments" / "{sample}_merged_aln_umi_annot_sorted_grouped.bam",
    output:
        bam=temp(wrkdir / "alignments" / "{sample}.cons.unmapped.bam"),
        metrics=logdir / "fgbio/call_consensus_reads.{sample}.log",
    params:
        min_reads=consensus_min_reads,
        min_base_qual=consensus_min_base_qual,
        min_input_base_mapq=consensus_min_input_base_mapq,
        error_rate_post_umi=consensus_error_rate_post_umi,
        error_rate_pre_umi=consensus_error_rate_pre_umi,
    threads: 4
    resources:
        mem_mb=4000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/fgbio.yaml"
    log:
        logdir / "fgbio/call_consensus_reads.{sample}.log",
    message:
        "Calling consensus reads from grouped reads."
    shell:
        "fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m --compression 0 CallMolecularConsensusReads "
        "--input {input.bam} "
        "--output {output.bam} "
        "--min-reads {params.min_reads} "
        "--error-rate-pre-umi {params.error_rate_pre_umi} "
        "--error-rate-post-umi {params.error_rate_post_umi} "
        "--min-consensus-base-quality {params.min_base_qual} "
        "--min-input-base-quality {params.min_input_base_mapq} "
        "--threads {threads} "
        "&> {log}"


rule sort_name_index:
    input:
        bam=wrkdir / "alignments" / "{sample}.cons.unmapped.bam",
    output:
        bam=temp(wrkdir / "alignments" / "{sample}.cons.unmapped.sorted.bam"),
    conda:
        "../envs/samtools.yaml"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        logdir / "samtools/{sample}_sort_name.log",
    message:
        "Sorting and indexing  concensus bam file"
    shell:
        " ( "
        " samtools sort --threads 8  -n -u -T {resources.tmpdir} -o {output.bam} {input.bam} "
        " ) &> {log} "


rule filter_consensus_reads:
    input:
        bam=wrkdir / "alignments" / "{sample}.cons.unmapped.sorted.bam",
        ref=genome,
    output:
        bam=temp(wrkdir / "alignments" / "{sample}.cons.filtered.bam"),
    params:
        min_reads=filter_min_reads,
        min_base_qual=filter_min_base_qual,
        max_base_error_rate=filter_max_base_error_rate,
        max_read_error_rate=filter_max_read_error_rate,
        max_no_call_fraction=filter_max_no_call_fraction,
    threads: 8
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/fgbio.yaml"
    log:
        logdir / "fgbio" / "filter_consensus_reads.{sample}.log",
    message:
        "Filtering consensus reads and sorting into coordinate order."
    shell:
        "("
        "fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m --compression 1 FilterConsensusReads "
        "--input {input.bam} "
        "--output {output.bam} "
        "--ref {input.ref} "
        "--min-reads {params.min_reads} "
        "--max-read-error-rate {params.max_read_error_rate} "
        "--max-base-error-rate {params.max_base_error_rate} "
        "--min-base-quality {params.min_base_qual} "
        "--max-no-call-fraction {params.max_no_call_fraction} "

        ") &> {log} "


rule gatherConsensusMetrics:
    input:
        logdir / "fgbio/call_consensus_reads.{sample}.log",
    output:
        wrkdir / "metrics" / "{sample}_consensus_metrics.tsv",
    params:
        sample=str(config["sample"]),
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        logdir / "metrics" / "{sample}_consensus_metrics.log",
    conda:
        "../envs/consensus.yaml"
    script:
        "../scripts/getConsensusMetrics.py"
