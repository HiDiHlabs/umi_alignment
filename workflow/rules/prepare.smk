
rule create_adapter_fastq:
    params:
        adapt_1=adapter_seq_r1,
        adapt_3=adapter_seq_r3
    output:
        adapt_1=temp(wrkdir / '{sample}' / 'cutadapt' / 'adapt_1.fastq'),
        adapt_3=temp(wrkdir / '{sample}' / 'cutadapt' / 'adapt_3.fastq'),
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

rule create_links_files:
    params:
        metadata=config['metadata'],
        fastq_dir=wrkdir / 'fastq'
    resources:
        mem_mb=1000,
        runtime=20,
        nodes=1
    output:
        fastq_r1=expand(wrkdir / 'fastq' / '{run_id}' / '{sample}_R1_{lane}.fastq.gz',filtered_product, run_id=RUN_ID, sample=config['sample'], lane=LANE),
        fastq_r2=expand(wrkdir / 'fastq' / '{run_id}' / '{sample}_R2_{lane}.fastq.gz',filtered_product, run_id=RUN_ID, sample=config['sample'], lane=LANE),
        fastq_r3=expand(wrkdir / 'fastq' / '{run_id}' / '{sample}_R3_{lane}.fastq.gz',filtered_product, run_id=RUN_ID, sample=config['sample'], lane=LANE)

    run:
        metadata = pd.read_csv(params.metadata)
        metadata = metadata[metadata['SAMPLE_NAME']==config['sample']]

        for index, row in metadata.iterrows():
            fastq_file=Path(row['FASTQ_FILE'])
            suffix='fastq'
            if '.gz' in row['FASTQ_FILE']:
                suffix+='.gz'
            output_file=params.fastq_dir / row['RUN_ID'] / (row['SAMPLE_NAME']+'_'+row['READ']+'_'+row['LANE_NO']+'.'+suffix) 
            os.makedirs(params.fastq_dir / row['RUN_ID'], exist_ok=True)
            if os.path.exists(output_file):
                os.remove(output_file)
            os.symlink(fastq_file, output_file)
                    
