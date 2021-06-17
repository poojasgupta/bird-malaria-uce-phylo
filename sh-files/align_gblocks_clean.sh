#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=align-gblocks-clean
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --constraint=EPYC             # node feature
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=250G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out       # Standard output log
#SBATCH --error=%x_%j.err        # Standard error log

#SBATCH --mail-user=pg84794@uga.edu
#SBATCH --mail-type=NONE

$SLURM_SUBMIT_DIR

module load phyluce/1.6.8	
	
	
time phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments /scratch/pg84794/UCE_run5/taxon-alignments/dataset2-mafft-nexus-internal-trimmed-gblocks \
    --output /scratch/pg84794/UCE_run5/taxon-alignments/dataset2-mafft-nexus-internal-trimmed-gblocks-clean \
    --cores 12
