rule fastqbam:
    """
    Converting fastq to bam to assign read group and library information
    This is important to ensure that RG and library information is added when alignment occurs
    becomes important for downstream analysis such as consensus calling
    """
    input:
        genome=genome,
        fastq_r1=(
            (
                wrkdir
                / "fastq"
                / "{run_id}"
                / "cutadapt"
                / "{sample}_R1_{lane}_trim.fastq.gz"
            )
            if config["trim_adapters"]
            else (wrkdir / "fastq" / "{run_id}" / "{sample}_R1_{lane}.fastq.gz")
        ),
        fastq_r3=(
            (
                (
                    wrkdir
                    / "fastq"
                    / "{run_id}"
                    / "cutadapt"
                    / "{sample}_R2_{lane}_trim.fastq.gz"
                )
                if config["trim_adapters"]
                else (wrkdir / "fastq" / "{run_id}" / "{sample}_R2_{lane}.fastq.gz")
            )
            if read_structure
            else (
                (
                    wrkdir
                    / "fastq"
                    / "{run_id}"
                    / "cutadapt"
                    / "{sample}_R3_{lane}_trim.fastq.gz"
                )
                if config["trim_adapters"]
                else (wrkdir / "fastq" / "{run_id}" / "{sample}_R3_{lane}.fastq.gz")
            )
        ),
    output:
        temp(wrkdir / "fastq" / "{run_id}" / "{sample}_{lane}_unmapped.bam"),
    params:
        library=library_prep_kit,
        read_structure="--read-structures " + read_structure if read_structure else "",
        read_group=lambda wc: (wc.run_id + "_" + wc.lane + "_" + wc.sample),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/fgbio.yaml"
    log:
        logdir / "fgbio" / "fastqtobam_{run_id}_{sample}_R1_{lane}.log",
    message:
        "Converting fastq to bam to assign read group and library information."
    shell:
        "("
        "fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m --compression 1 FastqToBam "
        "--input {input.fastq_r1} {input.fastq_r3} "
        "--sample {wildcards.sample} "
        "--library {params.library} "
        "--output {output} {params.read_structure} "
        "--read-group-id {params.read_group}"
        ") &> {log}"


if correct_umi:

    rule CorrectUMI:
        input:
            bam=wrkdir / "fastq" / "{run_id}" / "{sample}_{lane}_unmapped.bam",
            umi_file=umi_file,
        output:
            bam=temp(
                wrkdir
                / "fastq"
                / "{run_id}"
                / "{sample}_{lane}_unmapped_umi_corrected.bam"
            ),
            metrics=wrkdir
            / "metrics"
            / "correct_umi"
            / "{run_id}"
            / "{sample}_{lane}_umi_metrics.txt",
        params:
            max_mismatches=correct_umi_max_mismatches,
            min_distance=correct_umi_min_distance,
        threads: 1
        resources:
            mem_mb=2000,
            runtime=72 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        conda:
            "../envs/fgbio.yaml"
        log:
            logdir / "fgbio" / "correct_umi_{run_id}_{sample}_{lane}.log",
        message:
            "Correcting UMIs."
        shell:
            "("
            "fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m --compression 1 --async-io CorrectUmis "
            "--input {input.bam} "
            "--output {output.bam} "
            "--max-mismatches 3 "
            "--min-distance 2 "
            "--umi-files {input.umi_file} "
            "--metrics {output.metrics} "
            "--dont-store-original-umis ) &> {log}"


if not read_structure:
    print("Hi")

    rule AnnotateUMI:
        input:
            alignment=wrkdir / "alignments" / "{run_id}" / "{sample}_aln_{lane}.bam",
            fastq_r2=wrkdir / "fastq" / "{run_id}" / "{sample}_R2_{lane}.fastq.gz",
        output:
            temp(
                wrkdir
                / "alignments"
                / "{run_id}"
                / "{sample}_aln_{lane}_umi_annot.bam"
            ),
        threads: 1
        resources:
            mem_mb=8000,
            runtime=72 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        conda:
            "../envs/fgbio.yaml"
        log:
            logdir / "fgbio" / "annotate_umi_{run_id}_{sample}_{lane}.log",
        message:
            "Annotating BAM with UMIs from fastq file."
        shell:
            "fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m AnnotateBamWithUmis -i {input.alignment} -f {input.fastq_r2} -o {output} -t RX -q UQ -s true &> {log}"


rule bwa_map:
    """
    First pass alignemnt
    Aligning reads to the genome using BWA
    """
    input:
        genome=genome,
        bam=(
            wrkdir
            / "fastq"
            / "{run_id}"
            / "{sample}_{lane}_unmapped_umi_corrected.bam"
            if correct_umi
            else wrkdir / "fastq" / "{run_id}" / "{sample}_{lane}_unmapped.bam"
        ),
    output:
        temp(wrkdir / "alignments" / "{run_id}" / "{sample}_aln_{lane}.bam"),
    threads: 24
    resources:
        mem_mb=24000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/fgbio.yaml"
    message:
        "First pass alignemnt. Aligning reads to the genome using BWA."
    log:
        logdir / "bwa" / "first_pass_align_{run_id}_{sample}_{lane}.log",
    shell:
        "("
        "samtools fastq {input.bam} "
        "| bwa mem -Y -K 150000000 -t {threads} -p {input.genome} - "
        "| fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx4G --compression 1 --async-io ZipperBams "
        "--unmapped {input.bam} "
        "--ref {input.genome} "
        "--output {output} "
        ") &> {log}"


rule merge:
    """
    Merging bam files from different lanes/runs
    """
    input:
        expand(
            wrkdir / "alignments" / "{run_id}" / "{sample}_aln_{lane}.bam",
            filtered_product,
            run_id=RUN_ID,
            sample=config["sample"],
            lane=LANE,
        )
        if read_structure
        else expand(
            wrkdir / "alignments" / "{run_id}" / "{sample}_aln_{lane}_umi_annot.bam",
            filtered_product,
            run_id=RUN_ID,
            sample=config["sample"],
            lane=LANE,
        ),
    output:
        bam=temp(wrkdir / "alignments" / "{sample}_merged_umi_annot.bam"),
    threads: 24
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/samtools.yaml"
    message:
        "Merging bam files from different lanes/runs."
    log:
        logdir / "samtools/{sample}_merge.log",
    shell:
        "(samtools merge --threads {threads} -f {output.bam} {input}) &> {log} "


rule realign:
    """
    Second pass alignment using BWA once the consesnsus sequences called
    """
    input:
        bam=wrkdir / "alignments" / "{sample}.cons.filtered.bam",
        ref=genome,
    output:
        bam=temp(wrkdir / "alignments" / "{sample}.cons.filtered.realigned.bam"),
        bai=temp(wrkdir / "alignments" / "{sample}.cons.filtered.realigned.bam.bai"),
    threads: 28
    resources:
        mem_mb=80000,  # 8GB for BWA, 4GB for fgbio, 64GB for samtools sort and an overhead memory of 2GB
        runtime=72 * 60,
        nodes=1,
        mem_fgbio=4000,
        mem_samtools=8000,
        tmpdir=scratch_dir,
    params:
        samtools_threads=8,
        bwa_threads=24,
    conda:
        "../envs/fgbio.yaml"
    log:
        logdir / "bwa/{sample}_realign.log",
    message:
        "Second pass alignment using BWA on the consesnsus sequences called."
    shell:
        "("
        "samtools fastq {input.bam} "
        "| bwa mem -K 150000000 -Y -t {params.bwa_threads} -p {input.ref} - "
        "| fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_fgbio}m --compression 0 --async-io ZipperBams "
        "--unmapped {input.bam} "
        "--ref {input.ref} "
        "--tags-to-reverse Consensus "
        "--tags-to-revcomp Consensus "
        "| samtools sort --threads {params.samtools_threads} -m{resources.mem_samtools}m -T {resources.tmpdir} -o {output.bam}##idx##{output.bai} --write-index "
        ") &> {log} "
