# check if Mosdepth is run with in Exome/Panel or WGS mode

if seq_type in ["Panel", "WES"]:

    rule mosdepth:
        input:
            bam=wrkdir / "alignments" / "{sample}_dedup.recall.sorted.bam",
            target_regions=target_regions,
        output:
            out_1=wrkdir / "metrics" / "{sample}.mosdepth.global.dist.txt",
            out_2=wrkdir / "metrics" / "{sample}.mosdepth.summary.txt",
        log:
            logdir / "mosdepth/{sample}.log",
        conda:
            "../envs/mosdepth.yaml"
        threads: 1
        resources:
            mem_mb=8000,
            runtime=24 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        params:
            prefix=str(wrkdir / "metrics" / "{sample}"),
        message:
            "Running mosdepth for WES/panel data"
        shell:
            "mosdepth --by {input.target_regions} -n {params.prefix} {input.bam} &> {log}"

else:

    rule mosdepth:
        input:
            bam=wrkdir / "alignments" / "{sample}_dedup.recall.sorted.bam",
        output:
            out_1=wrkdir / "metrics" / "{sample}.mosdepth.global.dist.txt",
            out_2=wrkdir / "metrics" / "{sample}.mosdepth.summary.txt",
        log:
            logdir / "mosdepth/{sample}.log",
        conda:
            "../envs/mosdepth.yaml"
        threads: 1
        resources:
            mem_mb=8000,
            runtime=24 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        params:
            prefix=str(wrkdir / "metrics" / "{sample}"),
        message:
            "Running mosdepth for WGS data"
        shell:
            "mosdepth -n {params.prefix} {input.bam} &> {log}"
