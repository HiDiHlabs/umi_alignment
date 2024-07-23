# Deprecated


rule duplicates:
    input:
        wrkdir / "alignments" / "{sample}_temp.bam",
    output:
        bam=temp(wrkdir / "alignments" / "{sample}_dedup.bam"),
        # bai=temp(wrkdir / "alignments" / "{sample}_temp.bam.bai"),
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
        runtime=72 * 60,
        nodes=1,
    message:
        "Marking duplicates on pre-consensus reads.: Decrapated slated for removal."
    shell:
        " ( "
        " gatk --java-options '-Xmx{resources.mem_mb}m' MarkDuplicates -I {input} "
        " -O {output.bam} -M {output.metric} --BARCODE_TAG RX "
        " --SORTING_COLLECTION_SIZE_RATIO {params.SORTING_COLLECTION_SIZE_RATIO} "
        " ) &> {log} "
