print(
    """
    \tAlignment pipeline for UMI based sequencing reads
    \tAuthor: Shashwat Sahay
    \tEmail: shashwat.sahay@charite.de
    \tVersion: 0.2.0

    """
)


def get_data_time():
    now = datetime.now()

    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    return "========= " + dt_string + " ========= "


print(get_data_time(), "Starting run")

#######################################################################
################# Sanity check on the input data ######################
#######################################################################

if "sample" not in config:
    raise ValueError("Sample name not provided")

if "pid" not in config:
    raise ValueError("Patient ID not provided")

if "metadata" not in config:
    raise ValueError("Metadata file not provided")
else:
    metadata = pd.read_csv(config["metadata"])
    metadata = metadata[
        (metadata["SAMPLE_TYPE"] == config["sample"])
        & (metadata["PATIENT_ID"] == config["pid"])
    ]
    if metadata.shape[0] == 0:
        raise ValueError(
            "No metadata found for the given sample: %s and patient ID: %s"
            % (config["sample"], config["pid"])
        )

##########################################################
###### Setting Sequencing type and related parameters ####
##########################################################

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
            print(
                get_data_time(),
                "Bait regions not provided will be calculated from target regions",
            )
        else:
            bait_regions = config["bait_regions"]


##########################################################
################ Setting Library kit #####################
##########################################################

if "library_prep_kit" in config:
    library_prep_kit = config["library_prep_kit"]
    if library_prep_kit == "":
        library_prep_kit = "Unknown"
else:
    library_prep_kit = "Unknown"


##########################################################
################## Setting Genome ########################
##########################################################

if "genome" not in config:
    raise ValueError("Genome not provided")
else:
    genome = config["genome"]

if "dbsnp" not in config:
    raise ValueError("dbSNP file not provided")
else:
    dbsnp = config["dbsnp"]


##########################################################
############# Setting working directory ##################
##########################################################

if "work_dir" not in config:
    raise ValueError("Work directory not provided")
else:
    wrkdir = Path(config["work_dir"])

if "log_dir" not in config:
    logdir = wrkdir / "logs"
else:
    logdir = Path(config["log_dir"])

if "scratch_dir" not in config:
    scratch_dir = tempfile.gettempdir()
else:
    scratch_dir = config[
        "scratch_dir"
    ]  # snakemake doesnt accept path in resources has to be a string
    if not os.path.exists(scratch_dir):
        print(get_data_time(), "Scratch directory does not exist")
        os.makedirs(scratch_dir)

print(get_data_time(), "Setting working directory to %s" % wrkdir)
print(get_data_time(), "Setting log directory to %s" % logdir)
print(get_data_time(), "Setting temp directory to %s" % scratch_dir)


#########################################################################
############# Check if UMIs have been demultiplexed #####################
#########################################################################

if "read_structure" in config:
    print(
        get_data_time(),
        "Read structure has been provided expecting only R1 and R2 files for each sample and lane",
    )
    read_structure = config["read_structure"]
    print(read_structure)
    if read_structure == "":
        raise ValueError("Read structure defined but not provided")
else:
    read_structure = False


#########################################################################
############# Setting variables related to trimming #####################
#########################################################################

if "trim_adapters" not in config or (not config["trim_adapters"]):
    print(get_data_time(), "Adapter trimming has been turned off")
    config["trim_adapters"] = False
    adapter_seq_r1 = None
    adapter_seq_r3 = None
else:
    if config["trim_adapters"]:
        if "Adapter_R1" not in config:
            raise ValueError("Adapter sequence for R1 not provided")
        else:
            adapter_seq_r1 = config["Adapter_R1"]
        if "Adapter_R3" not in config:
            raise ValueError("Adapter sequence for R3 not provided")
        else:
            adapter_seq_r3 = config["Adapter_R3"]

##########################################################################
############# Setting variables related to UMI correction ################
##########################################################################

if "correct_umi" in config:
    correct_umi = config["correct_umi"]
else:
    correct_umi = False

if correct_umi:
    print(get_data_time(), "UMI Correction is enabled")
    if "umi_file" not in config:
        raise ValueError("UMI file not provided")
    else:
        umi_file = config["umi_file"]
    if "correct_umi_max_mismatches" in config:
        correct_umi_max_mismatches = config["correct_umi_max_mismatches"]
    else:
        print(
            get_data_time(), "Setting default value for correct_umi_max_mismatches to 3"
        )
        correct_umi_max_mismatches = 3
    if "correct_umi_min_distance" in config:
        correct_umi_min_distance = config["correct_umi_min_distance"]
    else:
        print(
            get_data_time(), "Setting default value for correct_umi_min_distance to 1"
        )
        correct_umi_min_distance = 1
else:
    print(get_data_time(), "UMI Correction is disabled")


##########################################################################
########### Initialising parameters for grouping by fgbio ################
##########################################################################

if "group_allowed_edits" not in config:
    print(get_data_time(), "Setting default value for allowed_edits to 1")
    group_allowed_edits = 1
else:
    group_allowed_edits = config["group_allowed_edits"]

if "group_strategy" not in config:
    print(get_data_time(), "Setting default value for strategy to Adjacency")
    group_strategy = "Adjacency"
else:
    group_strategy = config["group_strategy"]

if "group_min_mapq" not in config:
    print(get_data_time(), "Setting default value for group min mapq to 20")
    group_min_mapq = 20
else:
    group_min_mapq = config["group_min_mapq"]


##########################################################################
############# Intialising parameters for consensus calling ###############
##########################################################################

if "consensus_min_reads" not in config:
    print(get_data_time(), "Setting default value for consensus_min_reads to 1")
    consensus_min_reads = 1
else:
    consensus_min_reads = config["consensus_min_reads"]

if "consensus_min_base_qual" not in config:
    print(get_data_time(), "Setting default value for consensus_min_base_qual to 2")
    consensus_min_base_qual = 2
else:
    consensus_min_base_qual = config["consensus_min_base_qual"]


if "consensus_min_input_base_mapq" not in config:
    print(
        get_data_time(), "Setting default value for consensus_min_input_base_mapq to 10"
    )
    consensus_min_input_base_mapq = 10
else:
    consensus_min_input_base_mapq = config["consensus_min_input_base_mapq"]

if "consensus_error_rate_pre_umi" not in config:
    print(get_data_time(), "Setting default value for consensus_error_rate to 45")
    consensus_error_rate_pre_umi = 45
else:
    consensus_error_rate_pre_umi = config["consensus_error_rate_pre_umi"]

if "consensus_error_rate_post_umi" not in config:
    print(get_data_time(), "Setting default value for consensus_error_rate to 30")
    consensus_error_rate_post_umi = 30
else:
    consensus_error_rate_post_umi = config["consensus_error_rate_post_umi"]

############################################################################
######## Initialising parameters for consensus filtering ###################
############################################################################

if "filter_min_reads" not in config:
    print(get_data_time(), "Setting default value for filter_min_reads")
    filter_min_reads = 1
else:
    filter_min_reads = config["filter_min_reads"]

if "filter_min_base_qual" not in config:
    filter_min_base_qual = 10
else:
    filter_min_base_qual = config["filter_min_base_qual"]

print(get_data_time(), "Setting filter_min_base_qual to %s" % filter_min_base_qual)

if "filter_max_base_error_rate" not in config:
    print(
        get_data_time(), "Setting default value for filter_max_base_error_rate to 0.2"
    )
    filter_max_base_error_rate = 0.2
else:
    filter_max_base_error_rate = config["filter_max_base_error_rate"]

if "filter_max_read_error_rate" not in config:
    filter_max_read_error_rate = 0.1
    print(
        get_data_time(),
        "Setting default value for filter_max_base_error_rate to %s "
        % filter_max_read_error_rate,
    )
else:
    filter_max_read_error_rate = config["filter_max_read_error_rate"]

if "filter_max_no_call_fraction" not in config:
    print(
        get_data_time(), "Setting default value for filter_max_base_error_rate to 0.1"
    )
    filter_max_no_call_fraction = 0.1
else:
    filter_max_no_call_fraction = config["filter_max_no_call_fraction"]


# Create a set of valid combination for run lane and read ids

LANE = metadata["LANE_NO"].unique().tolist()
RUN_ID = metadata["RUN_ID"].unique().tolist()


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
