# UMI Based Sequencing Alignment Pipeline

Authour: Shashwat Sahay(shashwat.sahay@charite.de)

This pipeline was developed under supervision of Dr. Naveed Ishaque (naveed.ishaque@charite.de) and Prof. Roland Eils (roland.eils@charite.de).

The pipeline was tested and supported by Daniel Steiert.


The pipeline is made for aligning UMI based WGS/WES and Panel Seq data and to compute the QC metrics associated with it.

We require the sequencing is performed in paired end mode

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

# http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html
17. `group_allowed_edits`: Number of edit allowed when grouping based on umi. defaults to  0, should be set to zero if correct_umi is true
18. `group_min_mapq: 20`: Set `--min-map-q` of groupReadsByUMI
19. `group_strategy`: Set the `--strategy` param of groupReadsByUMI deafults to `Adjacency`

# https://fulcrumgenomics.github.io/fgbio/tools/latest/CallMolecularConsensusReads.html
20. `consensus_min_reads`: 1
21. `consensus_min_base_qual`: 2
22. `consensus_min_input_base_mapq`: 10
23. `consensus_error_rate_pre_umi`: 45
24. `consensus_error_rate_post_umi`: 30

# https://fulcrumgenomics.github.io/fgbio/tools/latest/FilterConsensusReads.html
25. `filter_min_reads`: 3
26. `filter_min_base_qual`: 2
27. `filter_max_base_error_rate`: 0.1
28. `filter_max_read_error_rate`: 0.05
29. `filter_max_no_call_fraction`: 0.2

# http://fulcrumgenomics.github.io/fgbio/tools/latest/FastqToBam.html
30. `read_structure`: 8M143T 8M143T

# http://fulcrumgenomics.github.io/fgbio/tools/latest/CorrectUmis.html
31. `correct_umi`: False
32. `correct_umi_max_mismatches`: 3
33. `correct_umi_min_distance`: 1
34. `umi_file`:
