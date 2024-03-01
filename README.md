# WGS_UMI_Alignment

The pipeline is made for aligning UMI based WGS data and to compute the QC metrics associated with it.

We require the sequencing is performed in paired end mode and must contain R1 (forward read) R2 (UMI) and R3 (reverse read) for each lane and run


# Setup/Installation

To run the pipeline make sure you have a working snakemake installation in a conda environment. We highly recommend using micromamba instead of any other alternatives conda!!!

Clone this repository with command
```
git clone https://github.com/HiDiHlabs/wgs_umi_alignment.git
```

And change to directory

```
cd wgs_umi_alignment
```


## Recommended Installation

Please create a conda environment

```
mamba env create -f workflow/envs/wgs-umi-dedup-base.yaml
```



## One step installation
### !!! Not Recommended !!!

A one for all conda environment is available at `workflow/envs/wgs-umi-dedup-full.yaml`. Although this is not the recommended way to prepare the conda environment in which the pipeline is run. As each rule has its own conda enviroment and can be/should be run independent of the base environment

```
mamba env create -f workflow/envs/wgs-umi-dedup-full.yaml
```


# Pipeline PREP
### !!! Important !!!

## Config file
To start the pipeline certain configurations must be made in the template config ```config/config.yaml```. It is recommended for each run of the pipeline a new config file be created based on the template

Please modify the entry for 
1. ```Adapter-Library-Prep-Kit_R1``` with the Adapter Sequences for Read 1 of the library prep. Needs to be `List`
2. ```Adapter-Library-Prep-Kit_R3``` with the Adapter Sequences for Read 1 of the library prep Needs to be `List`
3. ```metadata``` With the absolute path to the metadata sheet (Please check [Metadata section](#metadata) for format specifcation of the metadata sheet). Needs to be `Path`
4. ```work_dir``` Please provide with a work folder to store the output of pipeline. It is recommended that this folder is suffixed with PID_SAMPLE to avoid result overwriting. Needs to be `Path`
5. ```log_dir``` Please provide with a log folder to store the logs of the pipeline. It is recommended that this is inside the work_dir. Needs to be `Path`
6. ```sample``` Since this pipleine is run sample wise please mention the sample name as mentioned in the sample_name column of the metadata file. Needs to be `string`
7. ```pid``` The patient ID as that in the PATIENT ID column. Needs to be `string` 

## Metadata file

Please create a metadata file with columns
1. FASTQ_FILE: A column containing absolute paths to the locations of the FASTQ files,
2. READ: A column contain information with regards to R1, R2 and R3 of the sequencing file. Needs to be prefixed with `R` if not present
3. LANE_NO: A column containing the lane information for the sequencing files. Should be prefixed with `L_` if not present
4. SAMPLE_NAME: Containing the sample name which is inputed in the config file. Please note if the metadata file consists of multipe samples only the sample pid combination mentioned in the sample and pid directive of the config.yaml will be run
5. PATIENT_ID :Containing the ```pid``` which is inputed in the config file. Please note if the metadata file consists of multipe PIDs only the ```sample``` ```pid``` combination mentioned in the sample and pid directive of the config.yaml will be run.
6. LIB_PREP_KIT: Please inform the lib_prep_kit used for the sequencing 
7. RUN_ID: Please mention the run id for the squencing run for the sample 


# Running the pipeline 

Once you have the config file and the metadata file setup

First activate the conda environement containing base snakmake installation

```
mamba activate snakemake-wgs-umi-dedup-base
```

The run the pipeline with the following commands



```
 cd <Path/To/Pipeline/dir>
 snakemake --slurm -j 10 --configfile config/config.yaml --use-conda --conda-frontend mamba --profile profile 
```


### With one time install

#### !!! Not Recommended !!!
First activate the conda environement containing base snakmake installation

```
mamba activate snakemake-wgs-full
```

The run the pipeline with the following commands

``` 
cd <Path/To/Pipeline/dir>
snakemake --slurm -j 10 --configfile config/config.yaml --profile profile 
```

