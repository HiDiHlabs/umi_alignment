rule flagstatt_end:
    input:
        bam=os.path.join(wrkdir, alignment_dir, "{sample}_13-Sorted.bam"),
        bai=os.path.join(wrkdir, alignment_dir, "{sample}_13-Sorted.bam.bai"),
    output:
        os.path.join(wrkdir, metrics_dir, "{sample}_13-Sorted.flagstat"),
    conda:
        "../envs/sambamba.yaml"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        os.path.join(logdir, "sambamba/{sample}_13-Sorted.log"),
    message:
        "Running Flagstat"
    shell:
        "(sambamba flagstat -t {threads} {input.bam} > {output}) &> {log}"


rule flagstatt_primary_align:
    input:
        bam=os.path.join(wrkdir, alignment_dir, "{sample}_06-merged.bam"),
    output:
       os.path.join(wrkdir, metrics_dir, "{sample}_06-merged.flagstat"),
    conda:
        "../envs/sambamba.yaml"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        os.path.join(logdir, "sambamba/{sample}_06-merged.log"),
    message:
        "Running Flagstat"
    shell:
        "(sambamba flagstat -t {threads} {input.bam} > {output}) &> {log}"
