# disabling adapter trimming when when set to false,
# To maintain the same input and output files we create links instead of running cutadapt


rule create_adapter_fastq:
    params:
        adapt_1=adapter_seq_r1,
        adapt_2=adapter_seq_r2,
    output:
        adapt_1=temp(os.path.join(wrkdir, cutadapt_dir, "adapt_1.fastq")),
        adapt_2=temp(os.path.join(wrkdir, cutadapt_dir, "adapt_2.fastq")),
    threads: 1
    resources:
        mem_mb=1000,
        runtime=20,
        nodes=1,
        tmpdir=scratch_dir,
    message:
        "Creating adapter fastq files"
    run:
        with open(output.adapt_1, "w") as handle:
            count = 1
            for i in params.adapt_1:
                handle.write(">adapter_" + str(count) + "\n")
                handle.write(i + "\n")
                count += 1

        with open(output.adapt_2, "w") as handle:
            count = 1
            for i in params.adapt_2:
                handle.write(">adapter_" + str(count) + "\n")
                handle.write(i + "\n")
                count += 1

rule cutadapt:
    input:
        adapt_1=os.path.join(wrkdir, cutadapt_dir, "adapt_1.fastq"),
        adapt_2=os.path.join(wrkdir,  cutadapt_dir, "adapt_2.fastq"),
        fastq_r1=os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_R1_{lane}_00-softlink.part_{split}.fastq.gz"),
        fastq_r2=os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_R2_{lane}_00-softlink.part_{split}.fastq.gz"),
    output:
        fastq_r1=temp(os.path.join(wrkdir, fastq_dir , "{run_id}", cutadapt_dir, "{sample}_R1_{lane}_01-trim_{split}.fastq.gz")),
        fastq_r2=temp((os.path.join(wrkdir, fastq_dir, "{run_id}", cutadapt_dir, "{sample}_R2_{lane}_01-trim_{split}.fastq.gz"))),
        json_log = os.path.join(wrkdir, metrics_dir, "{run_id}_{sample}_{lane}_{split}_cutadapt_log.json"),
    params:
        cutadapt_params = cutadapt_params
    log:
        os.path.join(logdir, "cutadapt/{run_id}_{sample}_R1_R2_{lane}_{split}.log"),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/cutadapt.yaml"
    message:
        "Trimming adapters using cutadapt"
    shell:
        "cutadapt -j {threads} {params.cutadapt_params} --json={output.json_log} -a file:{input.adapt_1} -A file:{input.adapt_2} -o {output.fastq_r1} -p {output.fastq_r2} {input.fastq_r1} {input.fastq_r2} &> {log}"

rule subset_I1:
    input:
        fastq_r1=os.path.join(wrkdir, fastq_dir, "{run_id}", cutadapt_dir, "{sample}_R1_{lane}_01-trim_{split}.fastq.gz"),
        fastq_i1=os.path.join(wrkdir, fastq_dir, "{run_id}", "{sample}_I1_{lane}_00-softlink.fastq.gz"),
    output:
        fastq_i1=temp(os.path.join(wrkdir, fastq_dir, "{run_id}", "{sample}_I1_{lane}_00-softlink_{split}.fastq.gz")),
    log:
        os.path.join(logdir, "cutadapt/{run_id}_{sample}_R1_R2_{lane}_{split}.log"),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/pysam.yaml"
    message:
        "Subset I1 to match R1/R2"
    script:
        "../scripts/subset_I1.py"





