# check if Mosdepth is run with in Exome/Panel or WGS mode

if seq_type in ["Panel", "WES"]:

    rule mosdepth:
        input:
            bam=wrkdir / "alignments" / "{sample}_dedup.recall.sorted.bam",
            target_regions=target_regions,
        params:
            prefix=str(wrkdir / "metrics" / "{sample}"),
        output:
            out_1=wrkdir / "metrics" / "{sample}.mosdepth.global.dist.txt",
            out_2=wrkdir / "metrics" / "{sample}.mosdepth.summary.txt",
        threads: 1
        resources:
            mem_mb=8000,
            runtime=24 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        conda:
            "../envs/mosdepth.yaml"
        log:
            logdir / "mosdepth/{sample}.log",
        message:
            "Running mosdepth for WES/panel data"
        shell:
            "mosdepth --by {input.target_regions} -n {params.prefix} {input.bam} &> {log}"

else:

    rule mosdepth:
        input:
            bam=wrkdir / "alignments" / "{sample}_dedup.recall.sorted.bam",
        params:
            prefix=str(wrkdir / "metrics" / "{sample}"),
        output:
            out_1=wrkdir / "metrics" / "{sample}.mosdepth.global.dist.txt",
            out_2=wrkdir / "metrics" / "{sample}.mosdepth.summary.txt",
        threads: 1
        resources:
            mem_mb=8000,
            runtime=24 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        conda:
            "../envs/mosdepth.yaml"
        message:
            "Running mosdepth for WGS data"
        log:
            logdir / "mosdepth/{sample}.log",
        shell:
            "mosdepth -n {params.prefix} {input.bam} &> {log}"
