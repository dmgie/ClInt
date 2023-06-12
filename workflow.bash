## This is the main entry of the workflow bashscript file. We can call some sub-scripts here if 
## we want to make it a bit more modular.


## Set up various (environment) variables that we will use throughout the workflow
NUM_CORES=8
NUM_THREADS=8
MAX_RAM_GB=32

## Set up the paths to the various tools
SPADES_PATH=spades.py
BWA_PATH=bwa
SAMTOOLS_PATH=samtools
BCFTOOLS_PATH=bcftools
FREEBAYES_PATH=freebayes
VCFUTILS_PATH=vcfutils.pl
SNPEFF_PATH=snpeff.jar
SNPSIFT_PATH=SnpSift.jar
SNPEFF_CONFIG_PATH=snpEff.config
HUMAN_REFERENCE_PATH=human_g1k_v38.fasta
HUMAN_REFERENCE_GFF_PATH=human_g1k_v38.gff

############### .... a lot more needed

## Set up the paths to the various input files
READS_PATH=reads.fastq.gz










## Workflow
