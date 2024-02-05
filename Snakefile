from pathlib import Path
input_directory=Path('/dh-projects/richter_transformation/raw_data/WGS')
# input_directory=Path('/dh-projects/richter_transformation/raw_data/test_code')
output_directory=Path('/dh-projects/richter_transformation/analysis/results/WGS_alignment_cutadpt')
config['genome']="/applications/otp/reference-genomes/bwa06_1KGRef_PhiX/hs37d5_PhiX.fa",

adapter_seq_r1=['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'AGATCGGAAGAGCACACGTCTGAAC', 'TGGAATTCTCGGGTGCCAAGG', 'AGATCGGAAGAGCACACGTCT', 'CTGTCTCTTATACACATCT', 'AGATGTGTATAAGAGACAG']
adapter_seq_r3=['AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', 'AGATCGGAAGAGCGTCGTGTAGGGA', 'TGGAATTCTCGGGTGCCAAGG', 'AGATCGGAAGAGCACACGTCT', 'CTGTCTCTTATACACATCT', 'AGATGTGTATAAGAGACAG']



SAMPLES=['A006850326_211334_S1_L000','A006850326_211337_S2_L000','A006850326_211339_S3_L000','A006850326_211341_S4_L000','A006850326_211343_S5_L000','A006850326_211346_S6_L000']
# SAMPLES = ['test']
rule all:
    input:
        expand(output_directory / '{sample}' / 'metrics' / '{sample}.mosdepth.global.dist.txt',sample=SAMPLES),
        expand(output_directory / '{sample}' / "metrics" / 'flagstat' , sample=SAMPLES),
        expand(output_directory / '{sample}' / "metrics" / "insert_size_metrics.txt",sample=SAMPLES),
        expand(output_directory / '{sample}' / "metrics" / "insert_size.pdf",sample=SAMPLES),
        expand(output_directory / '{sample}' / 'metrics'/  '{sample}.mosdepth.summary.txt',sample=SAMPLES)

rule create_adapter_fastq:
    params:
        adapt_1=adapter_seq_r1,
        adapt_3=adapter_seq_r3
    output:
        adapt_1=output_directory / '{sample}' / 'cutadapt' / 'adapt_1.fastq',
        adapt_3=output_directory / '{sample}' / 'cutadapt' / 'adapt_3.fastq',
    threads: 1
    resources:
        mem_mb=1000,
        runtime=20,
        nodes=1
    run:
        with open(output.adapt_1, 'w') as handle:
            count=1
            for i in params.adapt_1:
                handle.write('>adapter_'+str(count)+'\n')
                handle.write(i+'\n')
                count+=1
                
        with open(output.adapt_3, 'w') as handle:
            count=1
            for i in params.adapt_3:
                handle.write('>adapter_'+str(count)+'\n')
                handle.write(i+'\n')
                count+=1
                
        
        
rule cutadapt:
    input:
        adapt_1=output_directory / '{sample}' / 'cutadapt' / 'adapt_1.fastq',
        adapt_3=output_directory / '{sample}' / 'cutadapt' / 'adapt_3.fastq',
        fastq_r1=input_directory / '{sample}_R1_001.fastq.gz',
        fastq_r3=input_directory / "{sample}_R3_001.fastq.gz"
    output:
        fastq_r1=output_directory / '{sample}' / 'cutadapt' / '{sample}_R1_001_trim.fastq.gz',
        fastq_r3=output_directory / "{sample}" / 'cutadapt' / '{sample}_R3_001_trim.fastq.gz',

    threads: 8
    resources:
        mem_mb=10000,
        runtime=4*60,
        nodes=1
    conda:
        "snakemake-richter"
    shell: "cutadapt -j 8 -a file:{input.adapt_1} -A file:{input.adapt_3} -o {output.fastq_r1} -p {output.fastq_r3} {input.fastq_r1} {input.fastq_r3}"



rule bwa_map:
    input:
        genome=config['genome'],
        fastq_r1=output_directory / '{sample}' / 'cutadapt' / '{sample}_R1_001_trim.fastq.gz',
        fastq_r3=output_directory / "{sample}" / 'cutadapt' / '{sample}_R3_001_trim.fastq.gz',
    output:
        output_directory / '{sample}' / "alignments" / "aln.bam"
    threads: 8
    resources:
        mem_mb=10000,
        runtime=24*60,
        nodes=1
    conda:
        "snakemake-richter"
    shell: "bwa mem -t {threads} -T 0 {input.genome} {input.fastq_r1} {input.fastq_r3} | samtools view -b -o {output}"

rule fgbio:
    input:
        alignment=output_directory / '{sample}' / "alignments" / "aln.bam",
        fastq_r2=input_directory / "{sample}_R2_001.fastq.gz"
    output:
        output_directory / '{sample}' / "alignments" / "umi_annot.bam"
    threads: 1
    resources:
        mem_mb=60000,
        runtime=4*60,
        nodes=1
    conda:
        "snakemake-richter"
    shell: "fgbio AnnotateBamWithUmis -i {input.alignment} -f {input.fastq_r2} -o {output} -t RX -q UQ -s true" 
    
rule sort_index:
    input:
        output_directory / '{sample}' / "alignments" /"umi_annot.bam"
    output:
        sort_bam=output_directory / '{sample}' / "alignments" / "umi_annot_sorted.bam",
        sort_bai=output_directory / '{sample}' / "alignments" /"umi_annot_sorted.bam.bai"
    threads: 8
    resources:
        mem_mb=10000,
        runtime=4*60,
        nodes=1
    conda:
        "snakemake-richter"
    shell: "samtools sort -o {output.sort_bam} -O bam -@ {threads} {input};samtools index -b {output.sort_bam} {output.sort_bai}"

rule duplicates:
    input:
        output_directory / '{sample}' / "alignments" / "umi_annot_sorted.bam"
    output:
        bam = output_directory / '{sample}' / "alignments" / "umi_annot_sorted_marked_duplicates.bam",
        bai = output_directory / '{sample}' / "alignments" / "umi_annot_sorted_marked_duplicates.bam.bai",
        metric = output_directory / '{sample}' / "metrics" / "marked_dup_metrics.txt"
    conda:
        "snakemake-richter"
    threads: 1
    resources:
        mem_mb=60000,
        runtime=24*60,
        nodes=1
    shell: "picard -Xmx60G MarkDuplicates -I {input} -O {output.bam} -M {output.metric} --BARCODE_TAG RX --SORTING_COLLECTION_SIZE_RATIO 0.1;samtools index -b {output.bam} {output.bai}"



rule flagstatt:
    input:
        bam=output_directory / '{sample}' / "alignments" / "umi_annot_sorted_marked_duplicates.bam"
    output:
        output_directory / '{sample}' / "metrics" / "flagstat"
    conda:
        "snakemake-richter"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=24*60,
        nodes=1
    shell: "sambamba flagstat {input.bam} > {output}"
    
        
rule mosdepth:
    input:
        bam=output_directory / '{sample}' / "alignments" / "umi_annot_sorted_marked_duplicates.bam"
    params:
        prefix=str(output_directory / '{sample}' / 'metrics'/ '{sample}' )
    output:
        out_1=output_directory / '{sample}' / 'metrics' / '{sample}.mosdepth.global.dist.txt',
        out_2=output_directory / '{sample}' / 'metrics'/  '{sample}.mosdepth.summary.txt'
    threads: 1
    resources:
        mem_mb=10000,
        runtime=24*60,
        nodes=1
    conda:
        "snakemake-richter"
        
    shell: "mosdepth -n {params.prefix} {input.bam}"

rule InsertSize:
    input:
        bam=output_directory / '{sample}' / "alignments" / "umi_annot_sorted_marked_duplicates.bam",
    output:
        size_metric=output_directory / '{sample}' / "metrics" / "insert_size_metrics.txt",
        pdf=output_directory / '{sample}' / "metrics" / "insert_size.pdf"
    conda:
        "snakemake-richter"
    threads: 1
    resources:
        mem_mb=60000,
        runtime=24*60,
        nodes=1
    shell: "picard -Xmx60G CollectInsertSizeMetrics -I {input.bam} -O {output.size_metric} -H {output.pdf}"