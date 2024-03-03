rule fastqbam:
    """
    Converting fastq to bam to assign read group and library information
    This is important to ensure that RG and library information is added when alignment occurs
    becomes important for downstream analysis such as consensus calling
    """
    input:
        genome=genome,
        fastq_r1=tmpdir / 'fastq' / '{run_id}' / 'cutadapt' / '{sample}_R1_{lane}_trim.fastq.gz',
        fastq_r3=tmpdir / 'fastq' / '{run_id}' / 'cutadapt' / '{sample}_R3_{lane}_trim.fastq.gz',
    output:
        temp(tmpdir / "fastq" / '{run_id}'/ '{sample}_{lane}_unmapped.bam')
    params:
        library=library_prep_kit
    threads: 8
    resources:
        mem_mb=8000,
        runtime=24*60,
        nodes=1
    conda:
        "../envs/fgbio.yaml"
    log:
        logdir / "fgbio" / "fastqtobam_{run_id}_{sample}_R1_{lane}.log"
    message: 
        "Converting fastq to bam to assign read group and library information."
    shell: 
        "("
        "fgbio -Xmx{resources.mem_mb}m --compression 1 FastqToBam "
        "--input {input.fastq_r1} {input.fastq_r3} "
        "--sample {wildcards.sample} "
        "--library {params.library} "
        "--output {output} "
        ") &> {log}"

rule bwa_map:
    """
    First pass alignemnt
    Aligning reads to the genome using BWA 
    """
    input:
        genome=genome,
        bam = tmpdir / "fastq" / '{run_id}'/ '{sample}_{lane}_unmapped.bam'
    output:
        temp(tmpdir / "alignments" / '{run_id}'/ '{sample}_aln_{lane}.bam')
    threads: 8
    resources:
        mem_mb=10000,
        runtime=24*60,
        nodes=1
    conda:
        "../envs/fgbio.yaml"
    message: 
        "First pass alignemnt. Aligning reads to the genome using BWA."
    log:
        logdir / "bwa" / "first_pass_align_{run_id}_{sample}_{lane}.log"
    shell: 
        "("
        "samtools fastq {input.bam} "
        "| bwa mem -Y -t {threads} -p {input.genome} - "
        "| fgbio -Xmx4G --compression 1 --async-io ZipperBams "
        "--unmapped {input.bam} "
        "--ref {input.genome} "
        "--output {output} "
        ") &> {log}"

rule merge:
    """
    Merging bam files from different lanes/runs
    """
    input:
        expand(tmpdir / "alignments" / '{run_id}'/ '{sample}_aln_{lane}_umi_annot.bam', filtered_product,  run_id=RUN_ID, sample=config['sample'], lane=LANE),
    output:
        temp(tmpdir / "alignments" / '{sample}_merged_umi_annot.bam'), 
    threads: 8
    resources:
        mem_mb=8000,
        runtime=24*60,
        nodes=1
    conda:
        "../envs/samtools.yaml"
    message: 
        "Merging bam files from different lanes/runs."
    log:
        logdir / "samtools/{sample}_merge.log"
    shell: "samtools merge -f {output} {input} &> {log}"


rule realign:
    """
    Second pass alignment using BWA once the consesnsus sequences called
    """
    input:
        bam = tmpdir / "alignments" / '{sample}.cons.filtered.bam',
        ref = genome,
    output:
        bam = temp(tmpdir / "alignments" / '{sample}.cons.filtered.realigned.bam'),
        bai = temp(tmpdir / "alignments" / '{sample}.cons.filtered.realigned.bam.bai'),
    threads: 8
    resources:
        mem_mb = 16000,
        runtime=24*60,
        nodes=1,
        mem_fgbio=8000,
    conda:
        "../envs/fgbio.yaml"
    log:
        logdir / "bwa/{sample}_realign.log"
    message: 
        "Second pass alignment using BWA on the consesnsus sequences called."
    shell: 
        "("
        "samtools fastq {input.bam} "
        "| bwa mem -Y -t {threads} -p {input.ref} - "
        "| fgbio -Xmx{resources.mem_fgbio}m --compression 0 --async-io ZipperBams "
        "--unmapped {input.bam} "
        "--ref {input.ref} "
        "--tags-to-reverse Consensus "
        "--tags-to-revcomp Consensus "
        "| samtools sort --threads 8 -o {output.bam}##idx##{output.bai} --write-index "
        ") &> {log} "
    
