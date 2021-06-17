#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=extract-fasta
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --constraint=EPYC             # node feature
#SBATCH --cpus-per-task=16             # Number of CPU cores per task
#SBATCH --mem=250G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out       # Standard output log
#SBATCH --error=%x_%j.err        # Standard error log

#SBATCH --mail-user=pg84794@uga.edu
#SBATCH --mail-type=NONE

cd $SLURM_SUBMIT_DIR

#create the data matrix configuration file

module load phyluce/1.6.8
#mkdir -p taxon-sets/dataset3

	
time phyluce_assembly_get_fastas_from_match_counts \
    --contigs /scratch/pg84794/UCE_run5/abyss-bam-en-assembly-k55/contigs/ \
    --locus-db /scratch/pg84794/UCE_run5/abyss-bam-en-uce-search-res/probe.matches.sqlite \
    --match-count-output /scratch/pg84794/UCE_run5/taxon-sets/dataset3/dataset3-taxa-incomp.conf \
    --output /scratch/pg84794/UCE_run5/taxon-sets/dataset3/dataset3-incomplete.fasta \
    --incomplete-matrix /scratch/pg84794/UCE_run5/taxon-sets/dataset3/dataset3-incomplete.incomplete \
    --extend-locus-db /scratch/pg84794/UCE_run5/metatranscript-uce-search-res/probe.matches.sqlite \
    --extend-locus-contigs /scratch/pg84794/UCE_run5/metatranscriptomes