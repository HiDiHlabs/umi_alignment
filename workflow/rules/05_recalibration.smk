rule baseRecalibrator:
    input:
        bam=os.path.join(wrkdir, alignment_dir, "{sample}_11-Realign-All-Reads.bam"),
        bai=os.path.join(wrkdir, alignment_dir, "{sample}_11-Realign-All-Reads.bam.bai"),
        dbsnp=dbsnp,
        genome=genome,
    output:
        table=os.path.join(wrkdir, metrics_dir, "{sample}_recal_data.table"),
    conda:
        "../envs/gatk.yaml"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        os.path.join(logdir, "gatk/{sample}_recal.log"),
    message:
        "Recalibrating with GATK BaseRecalibrator"
    shell:
        'gatk --java-options "-Djava.io.tmpdir={resources.tmpdir} -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_read_samtools=true -Xms4G -Xmx{resources.mem_mb}m -XX:ParallelGCThreads=2" BaseRecalibrator -I {input.bam} -R {input.genome} '
        ' --known-sites {input.dbsnp} '
        ' -O {output.table} &> {log} '


rule applyBSQR:
    input:
        bam=os.path.join(wrkdir, alignment_dir, "{sample}_11-Realign-All-Reads.bam"),
        bai=os.path.join(wrkdir, alignment_dir, "{sample}_11-Realign-All-Reads.bam.bai"),
        genome=genome,
        table=os.path.join(wrkdir, metrics_dir, "{sample}_recal_data.table"),
    output:
        bam=temp(os.path.join(wrkdir, alignment_dir, "{sample}_12-BaseRecalibrate.bam")),
        bai=temp(os.path.join(wrkdir, alignment_dir, "{sample}_12-BaseRecalibrate.bai")),
    conda:
        "../envs/gatk.yaml"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        os.path.join(logdir, "gatk/{sample}_applyBSQR.log"),
    message:
        "Recalibrating with GATK BaseRecalibrator"
    shell:
        'gatk --java-options "-Djava.io.tmpdir={resources.tmpdir} -Xms4G -Xmx{resources.mem_mb}m" '
        'ApplyBQSR --create-output-bam-index --emit-original-quals -I {input.bam} -R {genome} --bqsr-recal-file {input.table} -O {output.bam} &> {log} '


rule AnalyzeCovariates:
    input:
        table=os.path.join(wrkdir, metrics_dir, "{sample}_recal_data.table")
    output:
        analyse_covariates=os.path.join(wrkdir, metrics_dir, "{sample}_covariates.pdf"),
    conda:
        "../envs/gatk.yaml"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        os.path.join(logdir, "gatk/{sample}_analyzeBSQR.log"),
    message:
        "Recalibrating with GATK BaseRecalibrator"
    shell:
        "gatk AnalyzeCovariates "
        " -bqsr {input.table} "
        " -plots {output.analyse_covariates} "
        " &> {log} "


rule sort_final:
    input:
        bam=os.path.join(wrkdir, alignment_dir, "{sample}_12-BaseRecalibrate.bam"),
        bai=os.path.join(wrkdir, alignment_dir, "{sample}_12-BaseRecalibrate.bai"),
    output:
        bam=os.path.join(wrkdir, alignment_dir, "{sample}_13-Sorted.bam"),
    conda:
        "../envs/samtools.yaml"
    threads: 8
    params:
        mem_thread=8000,
    resources:
        mem_mb=8 * 8000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        os.path.join(logdir, "samtools/{sample}_sort.log"),
    message:
        "Sorting and indexing recalibrated bam file"
    shell:
        " samtools sort --threads 8 -m{params.mem_thread}m -o {output.bam} " ##idx##{output.bai} 
        "{input.bam} -T {resources.tmpdir}"
        " &> {log} "

rule sort_index:
    input:
        bam=os.path.join(wrkdir, alignment_dir, "{sample}_13-Sorted.bam"),
    output:
        bai=os.path.join(wrkdir, alignment_dir, "{sample}_13-Sorted.bam.bai"),
    conda:
        "../envs/samtools.yaml"
    threads: 2
    params:
        mem_thread=8000,
    resources:
        mem_mb=8000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    message:
        "Sorting and indexing recalibrated bam file"
    shell:
        " samtools index -b --threads {threads} -o {output.bai} " ##idx##{output.bai} 
        "{input.bam}"
