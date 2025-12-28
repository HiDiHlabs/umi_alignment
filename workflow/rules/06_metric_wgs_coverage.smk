rule coveragePlot:
    input:
        bam=wrkdir / "alignments" / "{sample}_13-Sorted.bam",
        bai=os.path.join(wrkdir, alignment_dir, "{sample}_13-Sorted.bam.bai"),
    params:
        binsize=50,
    output:
        plot=wrkdir / "metrics" / "{sample}_coverage.png",
    threads: 20
    resources:
        mem_mb=60000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/coveragePlot.yaml"
    log:
        logdir / "coveragePlot/{sample}.log",
    message:
        "Plotting Whole Genome coverage"
    script:
        "../scripts/getCoveragePlot_snakemake.R"
