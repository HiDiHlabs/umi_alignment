rule fastqbam:
    """
    Converting fastq to bam to assign read group and library information
    This is important to ensure that RG and library information is added when alignment occurs
    becomes important for downstream analysis such as consensus calling
    """
    input:
        genome=genome,
        fastq_r1=(os.path.join(wrkdir, fastq_dir , "{run_id}", cutadapt_dir, "{sample}_R1_{lane}_01-trim_{split}.fastq.gz") if config["trim_adapters"] 
                  else os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_R1_{lane}_00-softlink.part_{split}.fastq.gz")),
        fastq_r2=(os.path.join(wrkdir, fastq_dir , "{run_id}", cutadapt_dir, "{sample}_R2_{lane}_01-trim_{split}.fastq.gz") 
                  if config["trim_adapters"] else os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_R2_{lane}_00-softlink.part_{split}.fastq.gz")),
    output:
        temp(os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_{lane}_{split}_02-unmapped.bam") if not read_structure else os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_{lane}_{split}_03-unmapped_UMI-annot.bam")),
    params:
        library=library_prep_kit,
        read_structure="--read-structures " + read_structure if read_structure else "",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/fgbio.yaml"
    log:
        os.path.join(logdir, "fgbio", "fastqtobam_{run_id}_{sample}_R1_{lane}_{split}.log"),
    message:
        "Converting fastq to bam to assign read group and library information."
    shell:
        "("
        "fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m --compression 1 FastqToBam "
        "--input {input.fastq_r1} {input.fastq_r2} "
        "--sample {wildcards.sample} "
        "--library {params.library} "
        "--output {output} {params.read_structure} "
        ") &> {log}"


if not read_structure:
    rule AnnotateUMI:
        input:
            alignment=os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_{lane}_{split}_02-unmapped.bam"),
            fastq_i1=(os.path.join(wrkdir, fastq_dir, "{run_id}", "{sample}_I1_{lane}_00-softlink_{split}.fastq.gz")
                      if config["trim_adapters"] else os.path.join(wrkdir, fastq_dir, "{run_id}", "{sample}_I1_{lane}_00-softlink.fastq.gz")),
        output:
            temp(os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_{lane}_{split}_03-unmapped_UMI-annot.bam")),
        threads: 1
        resources:
            mem_mb=8000,
            runtime=72 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        conda:
            "../envs/fgbio.yaml"
        log:
            os.path.join(logdir, "fgbio", "annotate_umi_{run_id}_{sample}_{lane}_{split}.log"),
        message:
            "Annotating BAM with UMIs from fastq file."
        shell:
            "fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m AnnotateBamWithUmis "
            "-i {input.alignment} -f {input.fastq_i1} "
            "-o {output} -t RX -q UQ -s true --fail-fast true &> {log}" # -s true -> is this given? # Fail after first given UMI
            ## is not sorted anymore - as cutadapt might have removed some of the entries thus some UMIs are missing?!



if correct_umi:

    rule CorrectUMI:
        input:
            bam=os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_{lane}_{split}_03-unmapped_UMI-annot.bam"),
            umi_file=umi_file,
        output:
            bam=temp(os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_{lane}_{split}_03-unmapped_UMI-corrected.bam")),
            metrics=os.path.join(wrkdir, "metrics", "correct_umi", "{run_id}", "{sample}_{lane}_umi_metrics.txt"),
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
            os.path.join(logdir, "fgbio", "correct_umi_{run_id}_{sample}_{lane}_{split}.log"),
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


rule bwa_map:
    """
    First pass alignemnt
    Aligning reads to the genome using BWA
    """
    input:
        genome=genome,
        bam=(os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_{lane}_{split}_03-unmapped_UMI-corrected.bam") if correct_umi else os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_{lane}_{split}_03-unmapped_UMI-annot.bam")),
    output:
        bam = temp(os.path.join(wrkdir, alignment_dir, "{run_id}", "{sample}_{lane}_{split}_04-primary-aligned.bam")),
        # bai = temp(os.path.join(wrkdir, alignment_dir, "{run_id}", "{sample}_{lane}_{split}_04-primary-aligned.bam.bai")),
    threads: 25
    resources:
        mem_mb=12000,  # 8GB for BWA, 4GB for fgbio, 64GB for samtools sort and an overhead memory of 2GB
        runtime=72 * 60,
        nodes=1,
        mem_fgbio=4000,
        # mem_samtools=8000,
        tmpdir=scratch_dir,
    params:
        # samtools_threads=8,
        bwa_threads=24,
    conda:
        "../envs/fgbio.yaml"
    log:
        os.path.join(logdir, "bwa", "first_pass_align_{run_id}_{sample}_{lane}_{split}.log"),
    message:
        "First pass alignemnt. Aligning reads to the genome using BWA."
    shell:
        "("
        "samtools fastq {input.bam} "
        "| bwa mem -K 150000000 -Y -t {params.bwa_threads} -p {input.genome} - "
        "| fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_fgbio}m --compression 1 --async-io ZipperBams "
        "--unmapped {input.bam} "
        "--ref {input.genome} "
        "--output {output.bam} "
        # "| samtools sort --threads {params.samtools_threads} -m{resources.mem_samtools}m -T {resources.tmpdir} -o {output.bam}
        # "; samtools index --threads {threads} --bai --output {output.bai}  {output.bam} "
        ") &> {log}"



rule sortQueryName:
    """
    Downstream Tasks require sorting by QueryName 
    """
    input:
        bam=os.path.join(wrkdir, alignment_dir, "{run_id}", "{sample}_{lane}_{split}_04-primary-aligned.bam"),
    output:
        bam=temp(os.path.join(wrkdir, alignment_dir, "{run_id}", "{sample}_{lane}_{split}_05-QueryNameSorted.bam")),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/fgbio.yaml"
    log:
        os.path.join(logdir, "fgbio/querynameSort_{run_id}_{sample}_{lane}_{split}.log"),
    message:
        "Sort By QueryName"
    shell:
        "(fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m "
        "SortBam --input={input.bam} --sort-order=Queryname --output {output.bam}) >& log "


rule fix_mate:
    """
    Fixing mate information if required
    For some reason for some reads the mate information is not properly set.
    This can cause problems in downstream analysis.
    """
    input:
        bam=os.path.join(wrkdir, alignment_dir, "{run_id}", "{sample}_{lane}_{split}_05-QueryNameSorted.bam"),
    output:
        bam=temp(os.path.join(wrkdir, alignment_dir, "{run_id}", "{sample}_{lane}_{split}_05-mate-fix.bam")),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/fgbio.yaml"
    log:
        os.path.join(logdir, "fgbio/fixmate_{run_id}_{sample}_{lane}_{split}.log"),
    message:
        "Fixing mate information if required"
    shell:
        "(fgbio -Djava.io.tmpdir={resources.tmpdir} -Xmx{resources.mem_mb}m --compression 1 --async-io SetMateInformation "
        "--input {input.bam} "
        "--output {output.bam} "
        "--allow-missing-mates true )"
        " &> {log} "

rule merge:
    """
    Merging bam files from different lanes/runs
    """
    input:
        expand(os.path.join(wrkdir, alignment_dir, "{run_id}", "{sample}_{lane}_{split}_05-mate-fix.bam"), zip, **allow_dict_of_lists) #filtered_product, run_id=RUN_ID, sample=config["sample"], lane=LANE, split=split_list)
    output:
        bam=temp(os.path.join(wrkdir, alignment_dir, "{sample}_06-merged.bam")),
    threads: 6
    resources:
        mem_mb=16000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/samtools.yaml"
    message:
        "Merging bam files from different lanes/runs."
    log:
        os.path.join(logdir, "samtools/{sample}_merge.log"),
    shell:
        "(samtools merge --threads {threads} -f -o {output.bam} {input}) &> {log} "




rule realign:
    """
    Second pass alignment using BWA once the consesnsus sequences called
    """
    input:
        bam=os.path.join(wrkdir, alignment_dir, "{sample}_09-Consensus-Call-Filtered.bam"),
        ref=genome,
    output:
        bam=temp(os.path.join(wrkdir, alignment_dir, "{sample}_10-Realign-excl-supp.bam")),
        bai=temp(os.path.join(wrkdir, alignment_dir, "{sample}_10-Realign-excl-supp.bam.bai")),
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
        os.path.join(logdir, "bwa/{sample}_realign.log"),
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



rule extract_supp_reads:
    input:
        bam = os.path.join(wrkdir, alignment_dir, "{sample}_06-merged.bam"),
    output:
        bam = temp(os.path.join(wrkdir, alignment_dir, "{sample}_10-Supplemenary-Reads.bam")),
    threads: 2
    params:
        samtools_threads=1,
    resources:
        mem_mb=8000,
        mem_samtools=4000,
        nodes=1,
        runtime=1*24*60,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -f 0x800 -b {input.bam} | samtools sort --threads {params.samtools_threads} -m{resources.mem_samtools}m -o {output.bam}"



rule add_back_supp_reads:
    input:
        bam_supp = os.path.join(wrkdir, alignment_dir, "{sample}_10-Supplemenary-Reads.bam"),
        bam_realigned = os.path.join(wrkdir, alignment_dir, "{sample}_10-Realign-excl-supp.bam"),
    output:
        bam = temp(os.path.join(wrkdir, alignment_dir, "{sample}_11-Realign-All-Reads.bam")),
        bai = temp(os.path.join(wrkdir, alignment_dir, "{sample}_11-Realign-All-Reads.bam.bai")),
    threads: 1
    resources:
        mem_mb=4000,
        nodes=1,
        runtime=1*24*60,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools merge -o {output.bam} {input.bam_realigned} {input.bam_supp}; sleep 20s; samtools index -b {output.bam}"