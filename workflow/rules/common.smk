if "library_prep_kit" in config:
    library_prep_kit = config["library_prep_kit"]
    if library_prep_kit == "":
        library_prep_kit = "Unknown"
else:
    library_prep_kit = "Unknown"

if "sample" not in config:
    raise ValueError("Sample name not provided")

if "pid" not in config:
    raise ValueError("Patient ID not provided")

if "genome" not in config:
    raise ValueError("Genome not provided")
else:
    genome = config["genome"]

if "work_dir" not in config:
    raise ValueError("Work directory not provided")
else:
    wrkdir = Path(config["work_dir"])

if "log_dir" not in config:
    logdir = wrkdir / "logs"
else:
    logdir = Path(config["log_dir"])

if "trim_adapters" not in config:
    config["trim_adapters"] = False

if config["trim_adapters"]:
    if "Adapter_R1" not in config:
        raise ValueError("Adapter sequence for R1 not provided")
    else:
        adapter_seq_r1 = config["Adapter_R1"]
    if "Adapter_R3" not in config:
        raise ValueError("Adapter sequence for R3 not provided")
    else:
        adapter_seq_r3 = config["Adapter_R3"]

if "metadata" not in config:
    raise ValueError("Metadata file not provided")
else:
    metadata = pd.read_csv(config["metadata"])
    metadata = metadata[
        (metadata["SAMPLE_NAME"] == config["sample"])
        & (metadata["PATIENT_ID"] == config["pid"])
    ]
    if metadata.shape[0] == 0:
        raise ValueError("No metadata found for the given sample and patient ID")

if "SeqType" not in config:
    raise ValueError("Sequencing type not provided")
else:
    seq_type = config["SeqType"]
    if seq_type not in ["WGS", "WES", "Panel"]:
        raise ValueError("Invalid sequencing type provided")
    if seq_type in ["WES", "Panel"]:
        if "target_regions" not in config:
            raise ValueError("Target regions not provided")
        else:
            target_regions = config["target_regions"]

        if "chrom_sizes" not in config:
            raise ValueError("Chromosome sizes not provided")
        else:
            chrom_sizes = config["chrom_sizes"]
        if "dict_genome" not in config:
            raise ValueError("Genome not provided")
        else:
            dict_genome = config["dict_genome"]

        if "bait_regions" not in config:
            print("Bait regions not provided will be calculated from target regions")
        else:
            bait_regions = config["bait_regions"]

if "dbsnp" not in config:
    raise ValueError("dbSNP file not provided")
else:
    dbsnp = config["dbsnp"]

if "SORTING_COLLECTION_SIZE_RATIO" not in config:
    print("Setting default value for SORTING_COLLECTION_SIZE_RATIO to 0.01")
    SORTING_COLLECTION_SIZE_RATIO = 0.01
else:
    SORTING_COLLECTION_SIZE_RATIO = config["SORTING_COLLECTION_SIZE_RATIO"]

if "allowed_edits" not in config:
    print("Setting default value for allowed_edits to 1")
    allowed_edits = 1
else:
    allowed_edits = config["allowed_edits"]

if "strategy" not in config:
    print("Setting default value for strategy to Adjacency")
    strategy = "Adjacency"
else:
    strategy = config["strategy"]

if "consensus_min_reads" not in config:
    if seq_type == "WGS":
        print("Setting default value for consensus_min_reads to 3 as SeqType is WGS")
        consensus_min_reads = 3
    else:
        print(
            "Setting default value for consensus_min_reads to 1 as SeqType is not WGS"
        )
        consensus_min_reads = 1
else:
    consensus_min_reads = config["consensus_min_reads"]

if "consensus_min_base_qual" not in config:
    print("Setting default value for consensus_min_base_qual to 20")
    consensus_min_base_qual = 20
else:
    consensus_min_base_qual = config["consensus_min_base_qual"]

if "filter_min_reads" not in config:
    if seq_type == "WGS":
        print("Setting default value for filter_min_reads to 3 as SeqType is WGS")
        filter_min_reads = 3
    else:
        print("Setting default value for filter_min_reads to 1 as SeqType is not WGS")
        filter_min_reads = 1
else:
    filter_min_reads = config["filter_min_reads"]

if "filter_min_base_qual" not in config:
    print("Setting default value for filter_min_base_qual to 40")
    filter_min_base_qual = 40
else:
    filter_min_base_qual = config["filter_min_base_qual"]

if "filter_max_error_rate" not in config:
    print("Setting default value for filter_max_error_rate to 0.2")
    filter_max_error_rate = 0.2
else:
    filter_max_error_rate = config["filter_max_error_rate"]


LANE = metadata["LANE_NO"].unique().tolist()
RUN_ID = metadata["RUN_ID"].unique().tolist()


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
    for lane in metadata[metadata.RUN_ID == run_id]["LANE_NO"].unique().tolist():
        allow_list.add(
            frozenset(
                {("run_id", run_id), ("sample", config["sample"]), ("lane", lane)}
            )
        )


filtered_product = filter_combinator(product, allow_list)
