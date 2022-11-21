#!/bin/bash
#SBATCH -J Fastqc    # A single job name for the array
#SBATCH -n 4                       # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH --mem 35000                 # in MB
#SBATCH -t 0-06:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o job_%A_%a.log        # Standard output
#SBATCH -e job_%A_%a.log        # Standard error
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --account=mcgetthm-macrophage-fibroblast-rnaseq
#SBATCH --qos castles

set -e

module purge; module load bluebear

module load FastQC 0.11.9-Java-11


mkdir FASTQC

fastqc `find -maxdepth 5 -name '*fastq*' -print` --outdir=FASTQC