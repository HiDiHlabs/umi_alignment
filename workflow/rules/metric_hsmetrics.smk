# For some reason this conditional statement is needed to avoid a syntax error

if config['SeqType'] in ['Panel', 'WES']:
    rule createBait:
        input:
            capture_regions = config['capture_regions'],
            chrom_size = config['chrom_sizes'],
            dict_genome = config['dict_genome']
        output:
            flank_bed= temp(wrkdir / 'metrics' / '{sample}_flank.bed'),
            flank_intervals= temp(wrkdir / 'metrics' / '{sample}_flank.interval_list'),
            target_intervals= temp(wrkdir / 'metrics' / '{sample}_target.interval_list'),
        conda:
            "../envs/bedtools.yaml"
        threads: 1
        resources:
            mem_mb=8000,
            runtime=24*60,
            nodes=1
        log:
            logdir / "bedtools" / "{sample}_slop.log"
        message:
            "Creating flank and target intervals for Panel and WES data"
        shell: 
            "("
            " bedtools slop -b 100 -i {input.capture_regions} -g {input.chrom_size} > {output.flank_bed} "
            " && gatk BedToIntervalList -I {output.flank_bed} -O {output.flank_intervals} -SD {input.dict_genome}"
            " && gatk BedToIntervalList -I {input.capture_regions} -O {output.target_intervals} -SD {input.dict_genome}"
            ") &> {log}"

    rule HSmetrics:
        input:
            bam = wrkdir / "alignments" / '{sample}_dedup.recall.sorted.bam',
            bait_intervals = wrkdir / 'metrics' / '{sample}_flank.interval_list',
            target_intervals = wrkdir / 'metrics' / '{sample}_target.interval_list',
            genome=genome
        output:
            target=wrkdir / 'metrics' / '{sample}.hs_metrics.txt'
        conda:
            "../envs/gatk.yaml"
        threads: 1
        resources:
            mem_mb=8000,
            runtime=24*60,
            nodes=1
        message:
            "Calculating HS metrics for Panel and WES data"
        log:
            logdir / "picard/{sample}.hs_metrics.log"
        shell: "gatk CollectHsMetrics -I {input.bam} -O {output.target} -R {input.genome} -TI {input.target_intervals} -BI {input.bait_intervals} &> {log}"