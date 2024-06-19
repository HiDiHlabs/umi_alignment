# disabling adapter trimming when when set to false,
# To maintain the same input and output files we create links instead of running cutadapt

if config["trim_adapters"]:

    rule create_adapter_fastq:
        params:
            adapt_1=adapter_seq_r1,
            adapt_3=adapter_seq_r3,
        output:
            adapt_1=temp(wrkdir / "{sample}" / "cutadapt" / "adapt_1.fastq"),
            adapt_3=temp(wrkdir / "{sample}" / "cutadapt" / "adapt_3.fastq"),
        threads: 1
        resources:
            mem_mb=1000,
            runtime=20,
            nodes=1,
        message:
            "Creating adapter fastq files"
        run:
            with open(output.adapt_1, "w") as handle:
                count = 1
                for i in params.adapt_1:
                    handle.write(">adapter_" + str(count) + "\n")
                    handle.write(i + "\n")
                    count += 1

            with open(output.adapt_3, "w") as handle:
                count = 1
                for i in params.adapt_3:
                    handle.write(">adapter_" + str(count) + "\n")
                    handle.write(i + "\n")
                    count += 1

    rule cutadapt:
        input:
            adapt_1=wrkdir / "{sample}" / "cutadapt" / "adapt_1.fastq",
            adapt_3=wrkdir / "{sample}" / "cutadapt" / "adapt_3.fastq",
            fastq_r1=wrkdir / "fastq" / "{run_id}" / "{sample}_R1_{lane}.fastq.gz",
            fastq_r3=wrkdir / "fastq" / "{run_id}" / "{sample}_R3_{lane}.fastq.gz",
        output:
            fastq_r1=temp(
                wrkdir
                / "fastq"
                / "{run_id}"
                / "cutadapt"
                / "{sample}_R1_{lane}_trim.fastq.gz"
            ),
            fastq_r3=temp(
                wrkdir
                / "fastq"
                / "{run_id}"
                / "cutadapt"
                / "{sample}_R3_{lane}_trim.fastq.gz"
            ),
        log:
            logdir / "cutadapt/{run_id}_{sample}_R1_{lane}.log",
        threads: 8
        resources:
            mem_mb=8000,
            runtime=4 * 60,
            nodes=1,
        conda:
            "../envs/cutadapt.yaml"
        message:
            "Trimming adapters using cutadapt"
        shell:
            "cutadapt -j {threads} -a file:{input.adapt_1} -A file:{input.adapt_3} -o {output.fastq_r1} -p {output.fastq_r3} {input.fastq_r1} {input.fastq_r3} &> {log}"

else:

    rule cutadapt:
        input:
            fastq_r1=wrkdir / "fastq" / "{run_id}" / "{sample}_R1_{lane}.fastq.gz",
            fastq_r3=wrkdir / "fastq" / "{run_id}" / "{sample}_R3_{lane}.fastq.gz",
        output:
            fastq_r1=temp(
                wrkdir
                / "fastq"
                / "{run_id}"
                / "cutadapt"
                / "{sample}_R1_{lane}_trim.fastq.gz"
            ),
            fastq_r3=temp(
                wrkdir
                / "fastq"
                / "{run_id}"
                / "cutadapt"
                / "{sample}_R3_{lane}_trim.fastq.gz"
            ),
        log:
            logdir / "cutadapt/{run_id}_{sample}_R1_{lane}.log",
        threads: 1
        resources:
            mem_mb=1000,
            runtime=60,
            nodes=1,
        message:
            "Skipping adapter trimming, creating links instead of running cutadapt"
        conda:
            "../envs/cutadapt.yaml"
        shell:
            "(ln -s {input.fastq_r1} {output.fastq_r1} && ln -s {input.fastq_r3} {output.fastq_r3}) &> {log}"
