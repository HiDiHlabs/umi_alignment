
wrkdir = Path(config['work_dir'])
logdir = Path(config['log_dir'])

metadata = pd.read_csv(config['metadata'])
metadata = metadata[(metadata['SAMPLE_NAME']==config['sample']) & (metadata['PATIENT_ID']==config['pid'])]

lib_prep_kit = metadata['LIB_PREP_KIT'].tolist()[0]
if len(metadata['LIB_PREP_KIT'].unique())>1:
    raise ValueError('Multiple library prep kits detected. Please check yuor metadata file')

adapter_seq_r1 = config['Adapter-Library-Prep-Kit_R1']
adapter_seq_r3 = config['Adapter-Library-Prep-Kit_R3']


genome = config['genome'] 



LANE=metadata['LANE_NO'].unique().tolist()
RUN_ID=metadata['RUN_ID'].unique().tolist()


# Create a function which checks if the input configuration is valid

def filter_combinator(combinator, allow_list):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
            if frozenset(wc_comb) in allow_list:
                # print(wc_comb)
                yield wc_comb
    return filtered_combinator

allow_list = set()
for run_id in RUN_ID:
    for lane in metadata[metadata.RUN_ID==run_id]['LANE_NO'].unique().tolist():
        allow_list.add(frozenset({("run_id", run_id), ('sample', config['sample']), ("lane", lane)}))

if not config['trim_adapters']:
    print('Adapter trimming is disabled')

if config['SeqType'] in ['Panel', 'WES']:
    print("WES/ Panel data detected. Will run HS metrics")
    if not config['capture_regions']:
        raise ValueError('Capture regions not provided')





filtered_product = filter_combinator(product, allow_list)


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
    message:
        "Creating links to fastq files"
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


if config['trim_adapters']:
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
        message:
            "Creating adapter fastq files"
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