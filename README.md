# ClInt

## What is this pipeline

This pipeline was created as a way to integrate multiple methods and tools for analysis of clinical data - specifically in regards to analysing COVID19-related data. This includes mapping RNA-seq reads to a reference genome, carrying out variant calling (on RNA-seq data) and running circular RNA and micro-RNA prediction. It utilises Nextflow in order to carry out most of the process in parallel and on demand when the required input files are available, reducing the run time dramatically. 

Provided additionally is a Guix manifest file, which allows for running temporary containers via `guix shell` or also generating singularity or docker images via the `guix pack --format=docker` and `guix pack --format=squashfs` commands.

## Installation
To download the pipeline, run
```bash
git clone https://github.com/dmgie/ClInt.git
```
(please note the capital "I" in case errors arise during `git clone`)

Please also ensure that `Nextflow` has been installed on the system.

Afterwards, the pipeline can be run using `nextflow run ./main.nf <...>` where `<...>` are various parameters outlined in `Instructions and Usage` below.

### Prerequisite
The pipeline utilises the following programs
- *CircRNA*: DCC (circtools)
- *MiRNA*: miRanda, TargetScan
- *Mapping*: STAR, hisat
- *Variant Calling*: GATK
- *Others*: Nextflow, Singularity

*Nextflow* - This can be installed from their website: (here)[https://www.nextflow.io/]

*Singularity* - The pipeline utilises various containers in order to supply the necessary programs. These are handled within the scripts (and each process) via `Singularity` containers.

### Instructions and Usage

To use the pipeline, all that is needed is to provide the necessary files in the command line when calling the Nextflow command. An example is given below.

```sh
nextflow run main.nf --reference_file ../ref_genome.fna --input_dir ../reads_data/ --gff_file ../ref_genome.gff  --output_dir ./results -profile singularity
```

These parameters can also be given in the form a YAML/JSON formatted nextflow configuration file and provided with the `-params-file` parameter (more information in (Nextflow's Documentation)[https://www.nextflow.io/docs/latest/config.html]). 

The output of the pipeline will then be provided in the foldering given to the `--output_dir` parameter.

Instead of providing an input folder of reads, it is also possible to provide a samplesheet via the `--samplesheet` parameter. This is formatted as a CSV in 3 columns with "sample,fastq_1,fastq_2" as the headers, where `sample` is an identifier for a set of paired-end reads.

*Running inside Singularity*:
All that is required to run the pipeline, is to provide either a local path to the main "ClInt.sif" singularity container, or a link to a library which hosts the container. To do this, inside the `nextflow.config` file, modify the singularity profile - specifically the `container` variable inside the `process` block. Either provide a relative (`file://<relative_path>`), absolute path (`file:///<absolute_path>`) or a link. If given a link, it will download the containers into a temporary/cache directory, defined by the system's `SINGULARITY_CACHEDIR` environment variable.



## To Be Implemented
- [] Submit singularity image to a accessible repo to avoid need of building container
- [] The current implementation of the circRNA section of the pipeline is still incomplete. It has been heavily adapted from `nf-core/circRNA` in its current state. It is possible to run the circRNA detection and miRNA prediction steps separately using the output of the mapping step of the pipeline.
