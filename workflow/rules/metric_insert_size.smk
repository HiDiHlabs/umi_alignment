
rule InsertSize:
    input:
        bam=wrkdir / "alignments" / '{sample}_dedup.recall.sorted.bam',
    output:
        size_metric=wrkdir / "metrics" / "{sample}_insert_size_metrics.txt",
        pdf=wrkdir / "metrics" / "{sample}_insert_size.pdf"
    conda:
        "../envs/gatk.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=24*60,
        nodes=1
    log:
        logdir / "picard/{sample}_insert_size.log"
    message: 
        "Collecting insert size metrics"
    shell: "gatk CollectInsertSizeMetrics -I {input.bam} -O {output.size_metric} -H {output.pdf} &> {log}"
