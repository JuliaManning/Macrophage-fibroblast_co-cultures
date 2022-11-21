#!/bin/bash
#SBATCH -J BBDuk_trim_loop      # A single job name for the array
#SBATCH -n 4                       # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH --mem 8000                # in MB
#SBATCH -t 0-15:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o job_%A_%a.log        # Standard output
#SBATCH -e job_%A_%a.log        # Standard error
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --account=mcgetthm-macrophage-fibroblast-rnaseq
#SBATCH --qos castles

set -e
module purge; module load bluebear # this line is required

module load BBMap/38.87-GCC-8.3.0
#module load R/4.0.0-foss-2020a

# Remove adapter contamination, polyA read throuhgh, and low quality reads

for i in *.fastq; do bbduk.sh in=$i out=trimmed_clean_$i \
ref=truseq_rna.fa.gz literal=AAAAAAAAAAAAAAAAAA \
k=13 \
ktrim=r \
useshortkmers=t \
mink=5 \
qtrim=r \
trimq=10 \
minlength=20; done 

# in = chosen input file
# out = file name for out put 
# ref = file of adapter sequences (and in our case, polyA tail)
# k = kmer size to use
# ktrim = r, means right-trimming 
# mink allows use of shorter kmers at the end of reads
# qtrim = quality trim setting r trims right side only
# trimq = to Q10 using the Phred algorithrm, 


# hdist - hamming distance, allows for mismatch
# flt = force trimming, set how mny you'd wnat it to trim
# maq = quality filtering, discards reads below threshold
# ktrim = N, masks all 35-mers - to do with entropy 