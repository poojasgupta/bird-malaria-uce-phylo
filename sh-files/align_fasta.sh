#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=test-align-fasta
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

module load phyluce/1.6.8
	
time phyluce_align_seqcap_align \
    --fasta /scratch/pg84794/UCE_run5/gblocks-test/dataset3-incomplete.fasta \
    --output /scratch/pg84794/UCE_run5/gblocks-test/test-dataset3-mafft-nexus-internal-trimmed \
    --taxa 36 \
    --aligner mafft \
    --cores 16 \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --max-divergence 0.4
