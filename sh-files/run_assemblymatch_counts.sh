#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=raw-taxon-set
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

$SLURM_SUBMIT_DIR

#create the data matrix configuration file

module load phyluce/1.6.8
mkdir -p taxon-sets/datasetraw

time phyluce_assembly_get_match_counts \
    --locus-db /scratch/pg84794/UCE_run5/abyss-raw-uce-search-res/probe.matches.sqlite \
    --taxon-list-config /scratch/pg84794/UCE_run5/conf-files/taxon-set.conf \
    --taxon-group 'datasetraw' \
    --output taxon-sets/datasetraw/datasetraw-taxa-incomp.conf \
    --incomplete-matrix
