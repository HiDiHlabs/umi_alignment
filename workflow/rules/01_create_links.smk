rule create_links_files:
    params:
        metadata=config["metadata"],
        read_structure=read_structure,
    resources:
        mem_mb=1000,
        runtime=20,
        nodes=1,
    output:
        fastq_r1=expand(os.path.join(wrkdir, fastq_dir, "{run_id}", "{sample}_R1_{lane}_00-softlink.fastq.gz"),
                        # filtered_product,
                        run_id=RUN_ID, sample=config["sample"], lane=LANE),
        fastq_r2=expand(os.path.join(wrkdir, fastq_dir, "{run_id}", "{sample}_R2_{lane}_00-softlink.fastq.gz"),
                        # filtered_product,
                        run_id=RUN_ID, sample=config["sample"], lane=LANE),
        fastq_i1=expand(os.path.join(wrkdir, fastq_dir, "{run_id}", "{sample}_I1_{lane}_00-softlink.fastq.gz"),
                        # filtered_product,
                        run_id=RUN_ID, sample=config["sample"], lane=LANE) if not read_structure else [],
    # log:
    #     logdir, "{sample}_link-samples.log",
    message:
        "Creating links to fastq files"
    run:
        import time
        metadata = pd.read_csv(params.metadata)
        metadata = metadata[
            (metadata["SAMPLE_NAME"] == config["sample"])
            # & (metadata["PATIENT_ID"] == config["pid"])
        ]
        if params.read_structure:
            if metadata["READ"].nunique() != 2:
                raise ValueError("Read structure is provided but R1 R2 and I1 provided")
        else:
            if metadata["READ"].nunique() != 3:
                raise ValueError(
                    "Read structure is not provided but R1 R2 and I1 not provided"
                )
        for index, row in metadata.iterrows():
            fastq_file = Path(row["FASTQ_FILE"])
            suffix = "fastq"
            if row["FASTQ_FILE"].endswith("gz"): ## Only allow .gz files!
                suffix += ".gz"
            output_file = os.path.join(wrkdir, fastq_dir, row["RUN_ID"], 
                                       f'{row["SAMPLE_NAME"]}_{row["READ"]}_{row["LANE_NO"]}_00-softlink.{suffix}')
            os.makedirs(os.path.join(wrkdir, fastq_dir, row["RUN_ID"]) , exist_ok=True)
            if os.path.exists(output_file):
                os.remove(output_file)
            os.symlink(fastq_file, output_file)
        time.sleep(20)


if n_splits in [0, 1, "0", "1"]:
    rule simlink2:
        input:
            read1 = os.path.join(wrkdir, fastq_dir, "{run_id}", "{sample}_R1_{lane}_00-softlink.fastq.gz"),
            read2 = os.path.join(wrkdir, fastq_dir, "{run_id}", "{sample}_R2_{lane}_00-softlink.fastq.gz"),
        params:
            out_dir = os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir),
        output:
            read1 = temp(os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_R1_{lane}_00-softlink.part_000.fastq.gz")),
            read2 = temp(os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir, "{sample}_R2_{lane}_00-softlink.part_000.fastq.gz")),
        resources:
            mem_mb = 4000,
            runtime = 60 * 24 * 5,
            nodes = 1,
        threads: 1
        shell:
           "ln -s {input.read1} {output.read1};ln -s {input.read2} {output.read2}"
else:
    rule split_reads:
        input:
            read1 = os.path.join(wrkdir, fastq_dir, "{run_id}", "{sample}_R1_{lane}_00-softlink.fastq.gz"),
            read2 = os.path.join(wrkdir, fastq_dir, "{run_id}", "{sample}_R2_{lane}_00-softlink.fastq.gz"),
        params:
            n_splits = n_splits,
            out_dir = os.path.join(wrkdir, fastq_dir, "{run_id}", split_dir),
        output:
            read1 = temp(expand(os.path.join(wrkdir, fastq_dir, "{{run_id}}", split_dir, "{{sample}}_R1_{{lane}}_00-softlink.part_{split}.fastq.gz"), split=split_list)),
            read2 = temp(expand(os.path.join(wrkdir, fastq_dir, "{{run_id}}", split_dir, "{{sample}}_R2_{{lane}}_00-softlink.part_{split}.fastq.gz"), split=split_list)),
        conda: "../envs/seqkit.yaml"
        resources:
            mem_mb = 4000,
            runtime = 60 * 24 * 5,
            nodes = 1,
        threads: 1
        shell:
           "seqkit split2 -1 {input.read1} -2 {input.read2} --extension '.gz' -p {params.n_splits} -O {params.out_dir} " # --by-part-prefix '{sample}_R{read}_{lane}_00'

