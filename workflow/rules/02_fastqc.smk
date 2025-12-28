rule fastqc:
    input:
        fastq=os.path.join(wrkdir, fastq_dir, "{run_id}", "{sample}_R{run}_{lane}_00-softlink.fastq.gz"),
    params:
        outdir = os.path.join(wrkdir, metrics_dir, "{run_id}")
    output:
        html = os.path.join(wrkdir, metrics_dir, "{run_id}", "{sample}_R{run}_{lane}_00-softlink_fastqc.html"),
        zip = os.path.join(wrkdir, metrics_dir, "{run_id}", "{sample}_R{run}_{lane}_00-softlink_fastqc.zip"),
    conda:
        "../envs/fastqc.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        logdir / "fastqc/{run_id}_{sample}_R{run}_{lane}.log",
    message:
        "Running fastqc"
    shell:
        "mkdir -p {params.outdir};fastqc {input.fastq} -o {params.outdir} &> {log}"
