# ClInt

## What is it

This pipeline was created as a way to integrate multiple methods and tools to analyse clinical data - specifically regarding breast cancer and COVID-related data. It utilises Nextflow in order to carry out most of the process in parallel and on demand as soon as files become available, reducing the run time dramatically. 

Althought not yet implemented, utilising the SLURM executor for running it on HPC's also helps dramatically decrease the runtimes for some of the programs (especially during Assembly-based steps).

Provided additionally is a Dockerfile as well as a build script (`build_image.sh`) which can be used to create a Docker image, as well as a Singularity image - for usage depending on what is easier.

## Installation

### Prerequisite
<List of Programs>

<Nextflow>

<Singularity & Docker>

### Instructions

To use the pipeline, all that is needed is to give the necessary files in the command line when calling the Nextflow command. An example is given below.

```sh
nextflow run main.nf --reference_file ../ref_genome.fna --input_dir ../reads_data/ --gff_file ../ref_genome.gff  --output_dir ./
```

*Running inside Docker*:
All that is required to do to run it within the given Docker container is to add the `-profile docker` command to the end of the `nextflow run ...` command. 

To ensure that it runs, be sure to have built an image from the Dockerfile / `build_script.sh`.

*Running inside Singularity*:
Similarly to the Docker instruction, although instead of the `docker` profile, we use the `singularity` profile.

This also requires having the Singularity image built, and available in the same directory as the `main.nf` file

## Usage


## To Be Implemented
- [] SLURM configuration
- [] Strandedness and direction configuration
- [] Submit singularity image to a repo or something to download from - prevent needing local build everytime
