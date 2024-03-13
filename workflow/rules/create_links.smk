rule create_links_files:
    params:
        metadata=config["metadata"],
        fastq_dir=wrkdir / "fastq",
    resources:
        mem_mb=1000,
        runtime=20,
        nodes=1,
    output:
        fastq_r1=expand(
            wrkdir / "fastq" / "{run_id}" / "{sample}_R1_{lane}.fastq.gz",
            filtered_product,
            run_id=RUN_ID,
            sample=config["sample"],
            lane=LANE,
        ),
        fastq_r2=expand(
            wrkdir / "fastq" / "{run_id}" / "{sample}_R2_{lane}.fastq.gz",
            filtered_product,
            run_id=RUN_ID,
            sample=config["sample"],
            lane=LANE,
        ),
        fastq_r3=expand(
            wrkdir / "fastq" / "{run_id}" / "{sample}_R3_{lane}.fastq.gz",
            filtered_product,
            run_id=RUN_ID,
            sample=config["sample"],
            lane=LANE,
        ),
    message:
        "Creating links to fastq files"
    run:
        metadata = pd.read_csv(params.metadata)
        metadata = metadata[
            (metadata["SAMPLE_NAME"] == config["sample"])
            & (metadata["PATIENT_ID"] == config["pid"])
        ]

        for index, row in metadata.iterrows():
            fastq_file = Path(row["FASTQ_FILE"])
            suffix = "fastq"
            if ".gz" in row["FASTQ_FILE"]:
                suffix += ".gz"
            output_file = (
                params.fastq_dir
                / row["RUN_ID"]
                / (
                    row["SAMPLE_NAME"]
                    + "_"
                    + row["READ"]
                    + "_"
                    + row["LANE_NO"]
                    + "."
                    + suffix
                )
            )
            os.makedirs(params.fastq_dir / row["RUN_ID"], exist_ok=True)
            if os.path.exists(output_file):
                os.remove(output_file)
            os.symlink(fastq_file, output_file)
