#!/bin/bash
#SBATCH -J Alignment    # A single job name for the array
#SBATCH -n 4                       # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH --mem 35000                 # in MB
#SBATCH -t 0-10:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o job_%A_%a.log        # Standard output
#SBATCH -e job_%A_%a.log        # Standard error
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --account=mcgetthm-macrophage-fibroblast-rnaseq
#SBATCH --qos castles

set -e

module purge; module load bluebear

module load STAR/2.7.2b-GCC-8.3.0
module load SAMtools/1.9-foss-2018b

#STAR aligner
for i in *.fastq; do STAR --runThreadN ${SLURM_NTASKS} \
	 --genomeDir /castles/nr/projects/m/mcgetthm-macrophage-fibroblast-rnaseq/GenomeBuild/Human \
	 --readFilesIn $i \
	 --outFilterType BySJout \
	 --outFilterMultimapNmax 20 \
	 --alignSJoverhangMin 8 \
	 --alignSJDBoverhangMin 1 \
	 --outFilterMismatchNmax 999 \
	 --outFilterMismatchNoverLmax 0.1 \
	 --alignIntronMin 20 \
	 --alignIntronMax 1000000 \
	 --alignMatesGapMax 1000000 \
	 --limitBAMsortRAM 16000000000 \
	 --outSAMattributes NH HI NM MD \
	 --outSAMtype BAM SortedByCoordinate \
	 --outFileNamePrefix $i ; done

#Sam tools statistics and index 

for bamfile in *.sortedByCoord.out.bam ; do samtools index ${bamfile} ; done

for bamfile in *.sortedByCoord.out.bam ; do samtools flagstat ${bamfile}; done

#for i in *sortedByCoord.out.bam; do samtools index $iAligned.sortedByCoord.out.bam; done

# Encode standard options given below

# outFilterType By SJout - reduces the number of ”spurious” junctions
# outFilterMultimapNmax 20 - max number of multiple alignments per read
# alignSJoverhangMin 8 - min overhang for unannotated junctions
# alignSJDBoverhangMin 1 - min overhang for annotated junctions 
# outFilterMismatchNmax 999 - maximum number of mismatches per pair, large number switches off this filter
# outFilterMismatchNoverLmax 0.04 - max number of mismatches per pair relative to read length: for 2x100b
## max number of mismatches is 0.04*200=8 for the paired read (Lexogen suggested, 0.1)
# alignIntronMin 20 - minimum intron length
# alignIntronMax 1000000 - maximum intron length
# alignMatesGapMax 1000000 - maximum genomic distance between mates
# limitBAMsortRAM 16000000000 - default: 0 int>=0: maximum available RAM (bytes) for sorting BAM. 
## If =0, it will be set to the genome index size. 0 value can only be used with –genomeLoad NoSharedMemory option.
# outSAMattributes NH HI NM MD - defined in SAM format specs
# outSAMtype BAM SortedByCoordinate - BAM file (bianry equivalent of SAM file) can be sorted or not by Coordiantes


