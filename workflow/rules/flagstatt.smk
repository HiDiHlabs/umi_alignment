rule flagstatt:
    input:
        bam=wrkdir / "alignments" / "{sample}_dedup.recall.sorted.bam",
    output:
        wrkdir / "metrics" / "{sample}.flagstat",
    conda:
        "../envs/sambamba.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=24 * 60,
        nodes=1,
    log:
        logdir / "sambamba/{sample}.log",
    message:
        "Running Flagstat"
    shell:
        "(sambamba flagstat {input.bam} > {output}) &> {log}"
