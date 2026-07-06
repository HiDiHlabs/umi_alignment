rule coveragePlot:
    input:
        bam=wrkdir / "alignments" / "{sample}_dedup.recall.sorted.bam",
    output:
        plot=wrkdir / "metrics" / "{sample}_coverage.png",
    log:
        logdir / "coveragePlot/{sample}.log",
    conda:
        "../envs/coveragePlot.yaml"
    threads: 20
    resources:
        mem_mb=60000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    params:
        binsize=50,
    message:
        "Plotting Whole Genome coverage"
    script:
        "../scripts/getCoveragePlot_snakemake.R"
