
rule duplicates:
    input:
        wrkdir / "alignments" / "{sample}.cons.filtered.realigned.bam",
    output:
        bam=wrkdir / "alignments" / "{sample}_dedup.bam",
        bai=wrkdir / "alignments" / "{sample}_dedup.bam.bai",
        metric=wrkdir / "metrics" / "{sample}_marked_dup_metrics.txt",
    params:
        SORTING_COLLECTION_SIZE_RATIO=SORTING_COLLECTION_SIZE_RATIO,
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
        " -O {output.bam} -M {output.metric} "
        " --SORTING_COLLECTION_SIZE_RATIO {params.SORTING_COLLECTION_SIZE_RATIO} && "
        " samtools index -b {output.bam} {output.bai} "
        " ) &> {log} "
