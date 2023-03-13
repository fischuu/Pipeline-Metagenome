# Snakebite-Metagenomics

## Background and Features
A Snakemake - Metagenomics pipeline that covers the main aspects of shotgun
metagenomics.

The pipeline consists of different modules, namely

* QC
* Metagenome assembly
  * Full assembly, using all samples
  * Co-assembly subsets
* Gene prediction
* Read alignments
* Binning
* Taxonomic annotation

## Requirements
All required tools are defined as singularity / docker container that are
downloaded when needed. Further, the pre-configured resource allocations are
aimed for a slurm queueing system, other options are naturally possible, but
require then more user interaction. 

The pipeline works out-of-the-box, if the system has installed

* Docker / singularity /apptainer
* slurm

## Installation
The pipeline can be installed by cloning the corresponding GitHub repository

```
git clone https://github.com/fischuu/Snakebite-Metagenomics.git
```

## Preparation of the pipeline
There are three essential configuration files that needs to be adjusted

* `Snakebite-Metagenomics_config.yaml`
* `Snakebite-Metagenomics_server_config.yaml`
* `run_Snakebite-Metagenomics.sh`

The pipeline is shipped with an example configuration that can be adjusted to the
local needs.

In additon to that, a sample sheet needs to be created, an example for it is
also shipped with the pipeline

* `sampleSheet.tsv`

## Starting the pipeline
When all configuration files are filled, you can start the pipeline by running
the start script like this

```
bash run_Snakebite-Metagenomics.sh
```

## More details
For more details, please visit the project wiki page:

https://github.com/fischuu/Snakebite-Metagenomics/wiki

It contains details to the runtime paramters, explains how to populate the above
mentioned files and also contains convenience script to create the samplesheet
programmatically.

Further, the pipeline has several entry points so that it does not need to run
as a whole.