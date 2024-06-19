# For some reason this conditional statement is needed to avoid a syntax error

if seq_type in ["Panel", "WES"]:

    rule createBait:
        input:
            target_regions=target_regions,
            chrom_size=chrom_sizes,
        output:
            flank_bed=temp(wrkdir / "metrics" / "{sample}_flank.bed"),
        conda:
            "../envs/bedtools.yaml"
        threads: 1
        resources:
            mem_mb=1000,
            runtime=4 * 60,
            nodes=1,
        log:
            logdir / "bedtools" / "{sample}_slop.log",
        message:
            "Bait regions not provided creating flanks"
        shell:
            "("
            " bedtools slop -b 100 -i {input.target_regions} -g {input.chrom_size} > {output.flank_bed} "
            ") &> {log}"

    rule Intervals:
        input:
            target_regions=target_regions,
            bait_regions=(
                bait_regions
                if "bait_regions" in config
                else wrkdir / "metrics" / "{sample}_flank.bed"
            ),
            dict_genome=dict_genome,
        output:
            bait_intervals=temp(wrkdir / "metrics" / "{sample}_flank.interval_list"),
            target_intervals=temp(wrkdir / "metrics" / "{sample}_target.interval_list"),
        conda:
            "../envs/gatk.yaml"
        threads: 1
        resources:
            mem_mb=8000,
            runtime=4 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        log:
            logdir / "bedtools" / "{sample}_slop.log",
        message:
            "Creating flank and target intervals for Panel and WES data"
        shell:
            "("
            "gatk BedToIntervalList -I {input.bait_regions} -O {output.bait_intervals} -SD {input.dict_genome}"
            " && gatk BedToIntervalList -I {input.target_regions} -O {output.target_intervals} -SD {input.dict_genome}"
            ") &> {log}"

    rule HSmetrics:
        input:
            bam=wrkdir / "alignments" / "{sample}_dedup.recall.sorted.bam",
            bait_intervals=wrkdir / "metrics" / "{sample}_flank.interval_list",
            target_intervals=wrkdir / "metrics" / "{sample}_target.interval_list",
            genome=genome,
        output:
            target=wrkdir / "metrics" / "{sample}.hs_metrics.txt",
        conda:
            "../envs/gatk.yaml"
        threads: 1
        resources:
            mem_mb=8000,
            runtime=24 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        message:
            "Calculating HS metrics for Panel and WES data"
        log:
            logdir / "picard/{sample}.hs_metrics.log",
        shell:
            "gatk CollectHsMetrics -I {input.bam} -O {output.target} -R {input.genome} -TI {input.target_intervals} -BI {input.bait_intervals} &> {log}"
