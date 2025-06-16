# UMI Based Sequencing Alignment Pipeline

Authour: Shashwat Sahay(shashwat.sahay@charite.de)

This pipeline was developed under supervision of Dr. Naveed Ishaque (naveed.ishaque@charite.de) and Prof. Roland Eils (roland.eils@charite.de).

The pipeline was tested and supported by Daniel Steiert.


The pipeline is made for aligning UMI based WGS/WES and Panel Seq data and to compute the QC metrics associated with it.

We require the sequencing is performed in paired end mode


# Prerequisites

To run the pipeline make sure you have a working snakemake installation in a conda environment. We highly recommend using miniforge3 instead of any other alternatives conda!!!
Please follow this guide on how to install [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

Clone this repository with command

```
git clone https://github.com/HiDiHlabs/umi_alignment.git
```

And change to directory

```

cd umi_alignment
```


## Recommended Installation

Please create a conda environment

```
mamba env create -f workflow/envs/umi-dedup-base.yaml
```



## One step installation
### !!! Not Recommended !!!

A one for all conda environment is available at `workflow/envs/umi-dedup-full.yaml`. Although this is not the recommended way to prepare the conda environment in which the pipeline is run. As each rule has its own conda enviroment and can be/should be run independent of the base environment

```
mamba env create -f workflow/envs/umi-dedup-full.yaml
```


# Pipeline Preparation

### !!! Important !!!

## Config file
To start the pipeline certain configurations must be made in the template config ```config/config.yaml```. It is recommended for each run of the pipeline a new config file be created based on the template. It is also remcommended that the config file is stored in the output folder

Please modify the entry for

1. `SeqType`: Should be either `Panel`, `WGS` or `WES`

2. `library_prep_kit`: Library prep kit used for preparing the sample. If not available will be set to Unknown

3. `pid`: The patient ID as that in the PATIENT ID column. Needs to be `string`

4. `sample`: Since this pipleine is run sample wise please mention the sample name as mentioned in the sample_name column of the metadata file. Needs to be `string`
7. `metadata` Absolute path to the metadata sheet (Please check [Metadata section](#metadata) for format specifcation of the metadata sheet). Needs to be `Path`

8. `work_dir` Please provide an absolute path to a working folder to store the output of pipeline. It is recommended that this folder is suffixed with PID_SAMPLE to avoid result overwriting. Needs to be `Path`

9. `log_dir` Please provide an absolute path to a log folder to store the logs of the pipeline. It is recommended that this is inside the work_dir. Needs to be `Path`, if not provided or left empty the reverts to default `<workdir>/logs`

10. `genome` Please provide an abosolute path to the `genome.fa` file please note the genome should be indexed for use with bwa mem and indexes should be in the same folder as the genome

11. `dbsnp`: Please provide path to a vcf file used for recalibration by BaseRecalibrator

12. `trim_adapters`: A boolean to switch on and off the adapter trimming using cutadapt. It is highly recommended that adapter trimming be carried out but can switched off in rare cases

5. `Adapter_R1` with the Adapter Sequences for Read 1 of the library prep. Needs to be `List` can be an empty list if `trim_adapters` switch is set to False

6. `Adapter_R3` with the Adapter Sequences for Read 3 of the library prep Needs to be `List` can be an empty list if `trim_adapters` switch is set to False

13. `target_regions`: Absolute path to the target regions, must be set when `SeqType` is `WES` or `Panel`

14. `bait_regions`: Absolute path to the bait regions, if unset and `SeqType` is `WES` or `Panel`. A slop of 100bp on the `target_regions` is computed and used as bait regions

15. `chrom_sizes`: An absolute path to chromosomals length for the given genomes, ignored if `SeqType` is `WGS`

16. `dict_genome`: An absolute path to dict file for the given genomes, ignored if `SeqType` is `WGS`

17. `group_allowed_edits`: Number of edit allowed when grouping based on umi. defaults to  0, should be set to zero if correct_umi is true
18. `group_min_mapq: 20`: Set `--min-map-q` of groupReadsByUMI
19. `group_strategy`: Set the `--strategy` param of groupReadsByUMI deafults to `Adjacency`


20. `consensus_min_reads`: 1
21. `consensus_min_base_qual`: 2
22. `consensus_min_input_base_mapq`: 10
23. `consensus_error_rate_pre_umi`: 45
24. `consensus_error_rate_post_umi`: 30


25. `filter_min_reads`: 3
26. `filter_min_base_qual`: 2
27. `filter_max_base_error_rate`: 0.1
28. `filter_max_read_error_rate`: 0.05
29. `filter_max_no_call_fraction`: 0.2


30. `read_structure`: 8M143T 8M143T


31. `correct_umi`: True
32. `correct_umi_max_mismatches`: 3
33. `correct_umi_min_distance`: 1
34. `umi_file`:



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

    group_allowed_edits = 1
    group_strategy = "Adjacency"

    group_min_mapq = 20



##########################################################################
############# Intialising parameters for consensus calling ###############
##########################################################################

    consensus_min_reads = 1
    consensus_min_base_qual = 2
    consensus_min_input_base_mapq = 10
    consensus_error_rate_pre_umi = 45
    consensus_error_rate_post_umi = 30

############################################################################
############### Parameters for consensus filtering #########################
############################################################################

    filter_min_reads = 1


    filter_min_base_qual = 10

    filter_max_base_error_rate = 0.2


    filter_max_read_error_rate = 0.1


    filter_max_no_call_fraction = 0.1


    max_coverage = 10000

## Metadata file

Please create a metadata file with columns other columns can exist but not required
1. FASTQ_FILE: A column containing absolute paths to the locations of the FASTQ files,
2. READ: A column contain information with regards to R1, R2 and R3 of the sequencing file. Needs to be prefixed with `R` if not present
3. LANE_NO: A column containing the lane information for the sequencing files. Should be prefixed with `L_` if not present
4. SAMPLE_TYPE: Containing the sample name which is inputed in the config file. Please note if the metadata file consists of multipe samples only the sample pid combination mentioned in the sample and pid directive of the config.yaml will be run
5. PATIENT_ID :Containing the ```pid``` which is inputed in the config file. Please note if the metadata file consists of multipe PIDs only the ```sample``` ```pid``` combination mentioned in the sample and pid directive of the config.yaml will be run.
7. RUN_ID: Please mention the run id for the squencing run for the sample


# Running the pipeline

Once you have the config file and the metadata file setup

First activate the conda environement containing base snakmake installation

```
mamba activate umi-dedup-base
```

The run the pipeline with the following commands

```
cd <Path/to/work_dir>
snakemake -j 10 --configfile <Path/to/config.yaml> \
 --use-conda --conda-frontend mamba \
 --conda-prefix <Path/to/local/conda_envs> \
 --profile <Path/to/pipeline_dir/profile> \
 --snakefile <Path/to/pipeline_dir/workflow/Snakefile

```

The `--conda-prefix` will install the conda environment required at the particular location specified, this will helpful in maintaining a single version across runs. This is the recommended way.

**!!!Note!!!**
Please make sure before running mulitple instances of the pipeline, to run this command for one sample or for some test data so as to setup the environments.
You can ignore the --conda-prefix command but is highly recommned to use ut

## Running with Singularity

It is recommend to run the pipeline using a singularity container when working on High Performance Cluster. For this you would require to start the pipeline with


```
cd <Path/to/work_dir>
snakemake -j 10 --configfile <Path/to/config.yaml> \
 --use-singularity --singularity-prefix <Path/to/local/Singularity> \
 --singularity-args "-B /Path/to/data_folder1/:/Path/to/data_folder1,/Path/to/data_folder2/:/Path/to/data_folder2" \
 --use-conda --conda-frontend mamba \
 --conda-prefix <Path/to/local/conda_envs> \
 --profile <Path/to/pipeline_dir/profile> \
 --snakefile <Path/to/pipeline_dir/workflow/Snakefile
```

The `--singularity-prefix` will install the singularity environment required at the particular location specified, this will helpful in maintaining a single version across runs. This is the recommended way.

**!!!Note!!!**
Please make sure before running mulitple instances of the pipeline run this command for one sample or for some test data so as to setup the singularity image and the conda environment inside them.

The `--singularity-args` must be used to bind the folders/files required by the pipeline
(**Hint:** the files and folders inside the `config/config.yaml` file the location of the fastq files from the metadata file will need to be bound).

The current working directory is automatically bound by snakemake. Please look at documentation at [Apptainer](https://apptainer.org/docs/user/latest/introduction.html) and [Snakmake documentation](https://snakemake.readthedocs.io/en/v7.32.3/snakefiles/deployment.html#running-jobs-in-containers)

**Side Note**
If you are using a HPC please talk to your cluster admin and make sure all compute nodes have a running singularity version, For the first run/setup of the conda environment and singularity container requires the Node to have a working internet connection.



## With one time install

#### !!! Not Recommended !!!
First activate the conda environement containing base snakmake installation

```
mamba activate umi-wgs-full
```

The run the pipeline with the following commands

```
cd <Path/to/work_dir>
snakemake --slurm -j 10 --configfile <Path/to/config.yaml> --profile <Path/to/pipeline_dir/profile> --snakefile <Path/to/pipeline_dir/workflow/nakefile
```
