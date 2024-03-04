rule fix_mate:
    """
    Fixing mate information if required
    For some reason for some reads the mate information is not properly set.
    This can cause problems in downstream analysis.
    """
    input:
        bam=wrkdir / "alignments" / "{sample}_merged_umi_annot.bam",
    output:
        bam=temp(wrkdir / "alignments" / "{sample}_mate_fix.bam"),
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
        logdir / "fgbio/fixmate_{sample}.log",
    message:
        "Fixing mate information if required"
    shell:
        "fgbio -Xmx{resources.mem_mb}m --compression 1 --async-io SetMateInformation "
        "--input {input.bam} "
        "--output {output.bam} "
        "--allow-missing-mates true"
        "&> {log} "


rule duplicates:
    input:
        wrkdir / "alignments" / "{sample}.cons.filtered.realigned.bam",
    output:
        bam=wrkdir / "alignments" / "{sample}_dedup.bam",
        bai=wrkdir / "alignments" / "{sample}_dedup.bam.bai",
        metric=wrkdir / "metrics" / "{sample}_marked_dup_metrics.txt",
    conda:
        "../envs/gatk.yaml"
    threads: 8
    log:
        logdir / "gatk/{sample}_dedup.log",
    resources:
        mem_mb=8000,
        runtime=24 * 60,
        nodes=1,
    message:
        "Marking duplicates on consensus reads."
    shell:
        " ( "
        " gatk --java-options '-Xmx{resources.mem_mb}m' MarkDuplicates -I {input} "
        " -O {output.bam} -M {output.metric} --BARCODE_TAG RX "
        " --SORTING_COLLECTION_SIZE_RATIO 0.01 && "
        " samtools index -b {output.bam} {output.bai} "
        " ) &> {log} "
