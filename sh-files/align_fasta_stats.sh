#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=align-fasta-stats
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

cd $SLURM_SUBMIT_DIR

module load Miniconda3/4.7.10
source activate pgCondaEnv
#module load phyluce/1.6.8
#this script doesnt work	
		
time /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_align_get_align_summary_data  \
     --alignments /scratch/pg84794/UCE_run5/taxon-sets/taxon-alignments/dataset3-mafft-nexus-internal-trimmed \
     --output-stats /scratch/pg84794/UCE_run5/taxon-sets/taxon-alignments/dataset3-alignment-stats.csv \
     --show-taxon-counts \
     --cores 16
