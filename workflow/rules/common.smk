if "library_prep_kit" in config:
    library_prep_kit = config['library_prep_kit']
    if library_prep_kit=='':
        library_prep_kit = 'Unknown'
else:
    library_prep_kit = 'Unknown'

if "sample" not in config:
    raise ValueError('Sample name not provided')

if "pid" not in config:
    raise ValueError('Patient ID not provided')

if "genome" not in config:
    raise ValueError('Genome not provided')
else:
    genome = config['genome']

if "work_dir" not in config:
    raise ValueError('Work directory not provided')
else:
    wrkdir = Path(config['work_dir'])
    
if "tmp_dir" not in config:
    raise ValueError('Tempory Files directory not provided')
else:
    tmpdir = Path(config['tmp_dir'])

if "log_dir" not in config:
    logdir = wrkdir / 'logs'
else:
    logdir = Path(config['log_dir'])

if "trim_adapters" not in config:
    config['trim_adapters'] = False

if config['trim_adapters']:
    if "Adapter_R1" not in config:
        raise ValueError('Adapter sequence for R1 not provided')
    else:
        adapter_seq_r1 = config['Adapter_R1']
    if "Adapter_R3" not in config:
        raise ValueError('Adapter sequence for R3 not provided')
    else:
        adapter_seq_r3 = config['Adapter_R3']

if "metadata" not in config:
    raise ValueError('Metadata file not provided')
else:
    metadata = pd.read_csv(config['metadata'])
    metadata = metadata[
        (metadata['SAMPLE_NAME']==config['sample']) & 
        (metadata['PATIENT_ID']==config['pid'])
        ]
    if metadata.shape[0]==0:
        raise ValueError('No metadata found for the given sample and patient ID')

if 'SeqType' not in config:
    raise ValueError('Sequencing type not provided')
else:
    seq_type = config['SeqType']
    if seq_type not in ['WGS', 'WES', 'Panel']:
        raise ValueError('Invalid sequencing type provided')
    if seq_type in ['WES', 'Panel']:
        if "target_regions" not in config:
            raise ValueError('Target regions not provided')
        else:
            target_regions = config['target_regions']

        if "chrom_sizes" not in config:
            raise ValueError('Chromosome sizes not provided')
        else:
            chrom_sizes = config['chrom_sizes']
        if "dict_genome" not in config:
            raise ValueError('Genome not provided')
        else:
            dict_genome = config['dict_genome']

        if "bait_regions" not in config:
            print('Bait regions not provided will be calculated from target regions')
        else:
            bait_regions = config['bait_regions']

if "dbsnp" not in config:
    raise ValueError('dbSNP file not provided')
else:
    dbsnp = config['dbsnp']
        




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





filtered_product = filter_combinator(product, allow_list)


rule create_links_files:
    params:
        metadata=config['metadata'],
        fastq_dir=tmpdir / 'fastq'
    resources:
        mem_mb=1000,
        runtime=20,
        nodes=1
    output:
        fastq_r1=expand(tmpdir / 'fastq' / '{run_id}' / '{sample}_R1_{lane}.fastq.gz',filtered_product, run_id=RUN_ID, sample=config['sample'], lane=LANE),
        fastq_r2=expand(tmpdir / 'fastq' / '{run_id}' / '{sample}_R2_{lane}.fastq.gz',filtered_product, run_id=RUN_ID, sample=config['sample'], lane=LANE),
        fastq_r3=expand(tmpdir / 'fastq' / '{run_id}' / '{sample}_R3_{lane}.fastq.gz',filtered_product, run_id=RUN_ID, sample=config['sample'], lane=LANE)
    message:
        "Creating links to fastq files"
    run:
        metadata = pd.read_csv(params.metadata)
        metadata = metadata[
            (metadata['SAMPLE_NAME']==config['sample']) & 
            (metadata['PATIENT_ID']==config['pid'])
            ]

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
            adapt_1=temp(tmpdir / '{sample}' / 'cutadapt' / 'adapt_1.fastq'),
            adapt_3=temp(tmpdir / '{sample}' / 'cutadapt' / 'adapt_3.fastq'),
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
