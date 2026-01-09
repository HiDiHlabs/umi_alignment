rule coveragePlot:
    input:
        bam=wrkdir / "alignments" / "{sample}_13-Sorted.bam",
        bai=os.path.join(wrkdir, alignment_dir, "{sample}_13-Sorted.bam.bai"),
    params:
        binsize=100, # Only certain steps available 1000, 500, 100, 50, 30, 15, 10, 5, 1
    output:
        plot=wrkdir / "metrics" / "{sample}_coverage.png",
    threads: 5
    resources:
        mem_mb=20000,
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
