# UMI Based Sequencing Alignment Pipeline

Authour: Shashwat Sahay(shashwat.sahay@charite.de)

This pipeline was developed under supervision of Dr. Naveed Ishaque (naveed.ishaque@charite.de) and Prof. Roland Eils (roland.eils@charite.de).

The pipeline was tested and supported by Daniel Steiert.


The pipeline is made for aligning UMI based WGS/WES and Panel Seq data and to compute the QC metrics associated with it.

We require the sequencing is performed in paired end mode and must contain R1 (forward read) R2 (UMI) and R3 (reverse read) for each lane and run

Currently the ability to provide support is limited.


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



## Metadata file

Please create a metadata file with columns
1. FASTQ_FILE: A column containing absolute paths to the locations of the FASTQ files,
2. READ: A column contain information with regards to R1, R2 and R3 of the sequencing file. Needs to be prefixed with `R` if not present
3. LANE_NO: A column containing the lane information for the sequencing files. Should be prefixed with `L_` if not present
4. SAMPLE_NAME: Containing the sample name which is inputed in the config file. Please note if the metadata file consists of multipe samples only the sample pid combination mentioned in the sample and pid directive of the config.yaml will be run
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
