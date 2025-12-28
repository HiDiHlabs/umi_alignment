rule sortQueryName:
    """
    Downstream Tasks require sorting by QueryName
    """
    input:
        bam=os.path.join(wrkdir, alignment_dir, "{sample}_06-merged.bam"),
    output:
        bam=temp(os.path.join(wrkdir, alignment_dir, "{sample}_06-QueryNameSorted.bam")),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/fgbio.yaml"
    log:
        os.path.join(logdir, "fgbio/querynameSort_{sample}.log"),
    message:
        "Sort By QueryName"
    shell:
        "(fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m "
        "SortBam --input={input.bam} --sort-order=Queryname --output {output.bam}) >& log "
  


rule fix_mate:
    """
    Fixing mate information if required
    For some reason for some reads the mate information is not properly set.
    This can cause problems in downstream analysis.
    """
    input:
        bam=os.path.join(wrkdir, alignment_dir, "{sample}_06-QueryNameSorted.bam"),
    output:
        bam=temp(os.path.join(wrkdir, alignment_dir, "{sample}_07-mate-fix.bam")),
    threads: 1
    resources:
        mem_mb=50000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/fgbio.yaml"
    log:
        os.path.join(logdir, "fgbio/fixmate_{sample}.log"),
    message:
        "Fixing mate information if required"
    shell:
        "(fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m --compression 1 --async-io SetMateInformation "
        "--input {input.bam} "
        "--output {output.bam} "
        "--allow-missing-mates true )"
        " &> {log} "


rule group_reads:
    input:
        bam=os.path.join(wrkdir, alignment_dir, "{sample}_07-mate-fix.bam"),
    output:
        bam=temp(os.path.join(wrkdir, alignment_dir, "{sample}_08-group-reads.bam")),
        stats=os.path.join(wrkdir, metrics_dir, "{sample}.grouped-family-sizes.txt"),
    params:
        strategy=group_strategy,
        allowed_edits=group_allowed_edits,
        min_mapq=group_min_mapq,
        # include_supplementary = "--include-supplementary" if include_supplementary_GroupReadsByUmi else "",
        # include_secondary = "--include-secondary" if include_secondary_GroupReadsByUmi else "",
    threads: 2
    resources:
        mem_mb=50000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/fgbio.yaml"
    log:
        os.path.join(logdir, "fgbio/group_{sample}.log"),
    message:
        "Grouping reads by UMI and position for consensus calling."
    shell:
        "fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m --compression 1 --async-io GroupReadsByUmi "
        "--input {input.bam} "
        "--strategy {params.strategy} "
        "--edits {params.allowed_edits} "
        "--min-map-q {params.min_mapq} "
        "--output {output.bam} "
        # "{params.include_secondary} {params.include_supplementary} "
        "--family-size-histogram {output.stats} "
        "&> {log} "

if not duplex_read_structure:
    rule call_filter_consensus_reads:
        input:
            bam=os.path.join(wrkdir, alignment_dir, "{sample}_08-group-reads.bam"),
            ref=genome,
        output:
            bam=temp(os.path.join(wrkdir, alignment_dir, "{sample}_09-Consensus-Call-Filtered.bam")),
            metrics=os.path.join(logdir, "fgbio/call_consensus_reads.{sample}.log"),
        params:
            min_reads=consensus_min_reads,
            min_input_base_mapq=consensus_min_input_base_mapq,
            error_rate_post_umi=consensus_error_rate_post_umi,
            error_rate_pre_umi=consensus_error_rate_pre_umi,
            filter_min_reads=filter_min_reads,
            min_base_qual=filter_min_base_qual,
            max_base_error_rate=filter_max_base_error_rate,
            max_read_error_rate=filter_max_read_error_rate,
            max_no_call_fraction=filter_max_no_call_fraction,
            memory_consensus=25000,
            memory_filter=25000,
        threads: 24
        resources:
            mem_mb=50000,
            runtime=72 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        conda:
            "../envs/fgbio.yaml"
        log:
            os.path.join(logdir, "fgbio/call_consensus_reads.{sample}.log"),
        message:
            "Calling consensus reads from grouped reads."
        shell:
            "( fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{params.memory_consensus}m --compression 0 CallMolecularConsensusReads "
            "--input {input.bam} "
            "--output /dev/stdout "
            "--min-reads {params.min_reads} "
            "--error-rate-pre-umi {params.error_rate_pre_umi} "
            "--error-rate-post-umi {params.error_rate_post_umi} "
            "--min-input-base-quality {params.min_input_base_mapq} "
            "--threads {threads} | fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{params.memory_filter}m --compression 1 FilterConsensusReads "
            "--input /dev/stdin "
            "--output {output.bam} "
            "--ref {input.ref} "
            "--min-reads {params.filter_min_reads} "
            "--max-read-error-rate {params.max_read_error_rate} "
            "--max-base-error-rate {params.max_base_error_rate} "
            "--min-base-quality {params.min_base_qual} "
            "--max-no-calls {params.max_no_call_fraction} ) "
            "&> {log}"

if duplex_read_structure:
    rule call_filter_duplex_consensus_reads:
        input:
            bam=os.path.join(wrkdir, alignment_dir, "{sample}_08-group-reads.bam"),
            ref=genome,
        output:
            bam=temp(os.path.join(wrkdir, alignment_dir, "{sample}_09-Consensus-Call-Filtered.bam")),
            metrics=os.path.join(logdir, "fgbio/call_consensus_reads.{sample}.log"),
        params:
            min_reads=consensus_min_reads,
            min_input_base_mapq=consensus_min_input_base_mapq,
            error_rate_post_umi=consensus_error_rate_post_umi,
            error_rate_pre_umi=consensus_error_rate_pre_umi,
            filter_min_reads=filter_min_reads,
            min_base_qual=filter_min_base_qual,
            max_base_error_rate=filter_max_base_error_rate,
            max_read_error_rate=filter_max_read_error_rate,
            max_no_call_fraction=filter_max_no_call_fraction,
            memory_consensus=25000,
            memory_filter=25000,
        threads: 24
        resources:
            mem_mb=50000,
            runtime=72 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        conda:
            "../envs/fgbio.yaml"
        log:
            os.path.join(logdir, "fgbio/call_consensus_reads.{sample}.log"),
        message:
            "Calling consensus reads from grouped reads."
        shell:
            "( fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{params.memory_consensus}m --compression 0 CallDuplexConsensusReads "
            "--input {input.bam} "
            "--output /dev/stdout "
            "--min-reads {params.min_reads} "
            "--error-rate-pre-umi {params.error_rate_pre_umi} "
            "--error-rate-post-umi {params.error_rate_post_umi} "
            "--min-input-base-quality {params.min_input_base_mapq} "
            "--threads {threads} | fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{params.memory_filter}m --compression 1 FilterConsensusReads "
            "--input /dev/stdin "
            "--output {output.bam} "
            "--ref {input.ref} "
            "--min-reads {params.filter_min_reads} "
            "--max-read-error-rate {params.max_read_error_rate} "
            "--max-base-error-rate {params.max_base_error_rate} "
            "--min-base-quality {params.min_base_qual} "
            "--max-no-call-fraction {params.max_no_call_fraction} ) "
            "&> {log}"

# rule sort_name_index:
#     input:
#         bam=os.path.join(wrkdir, alignment_dir, "{sample}.cons.unmapped.bam",
#     output:
#         bam=temp(os.path.join(wrkdir, alignment_dir, "{sample}.cons.unmapped.sorted.bam"),
#     conda:
#         "../envs/samtools.yaml"
#     threads: 8
#     params:
#         mem_thread=8000,
#     resources:
#         mem_mb=8 * 8000,
#         runtime=24 * 60,
#         nodes=1,
#         tmpdir=scratch_dir,
#     log:
#         logdir, "samtools/{sample}_sort_name.log",
#     message:
#         "Sorting and indexing  concensus bam file"
#     shell:
#         " ( "
#         " samtools sort --threads 8 -m{params.mem_thread}m -n -u -T {resources.tmpdir} -o {output.bam} {input.bam} "
#         " ) &> {log} "


# rule filter_consensus_reads:
#     input:
#         bam=os.path.join(wrkdir, alignment_dir, "{sample}.cons.unmapped#.bam",
#         ref=genome,
#     output:
#         bam=temp(os.path.join(wrkdir, alignment_dir, "{sample}.cons.filtered#.bam"),
#     params:
#         min_reads=filter_min_reads,
#         min_base_qual=filter_min_base_qual,
#         max_base_error_rate=filter_max_base_error_rate,
#         max_read_error_rate=filter_max_read_error_rate,
#         max_no_call_fraction=filter_max_no_call_fraction,
#     threads: 8
#     resources:
#         mem_mb=50000,
#         runtime=72 * 60,
#         nodes=1,
#         tmpdir=scratch_dir,
#     conda:
#         "../envs/fgbio.yaml"
#     log:
#         logdir, "fgbio", "filter_consensus_reads.{sample}.log",
#     message:
#         "Filtering consensus reads and sorting into coordinate order."
#     shell:
#         "("
#         "fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m --compression 1 FilterConsensusReads "
#         "--input {input.bam} "
#         "--output {output.bam} "
#         "--ref {input.ref} "
#         "--min-reads {params.min_reads} "
#         "--max-read-error-rate {params.max_read_error_rate} "
#         "--max-base-error-rate {params.max_base_error_rate} "
#         "--min-base-quality {params.min_base_qual} "
#         "--max-no-call-fraction {params.max_no_call_fraction} "

#         ") &> {log} "


rule gatherConsensusMetrics:
    input:
        os.path.join(logdir, "fgbio/call_consensus_reads.{sample}.log"),
    output:
        os.path.join(wrkdir, metrics_dir, "{sample}_consensus_metrics.tsv"),
    params:
        sample=str(config["sample"]),
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        os.path.join(logdir, metrics_dir, "{sample}_consensus_metrics.log"),
    conda:
        "../envs/consensus.yaml"
    script:
        "../scripts/getConsensusMetrics.py"
