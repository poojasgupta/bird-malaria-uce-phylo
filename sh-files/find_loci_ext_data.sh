#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=ext-data-genomes-transcriptomes-editnames-uce-search65
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --constraint=EPYC             # node feature
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=400G
#SBATCH --time=100:00:00
#SBATCH --output=%x_%j.out       # Standard output log
#SBATCH --error=%x_%j.err        # Standard error log

#SBATCH --mail-user=pg84794@uga.edu
#SBATCH --mail-type=NONE

cd $SLURM_SUBMIT_DIR

module load Miniconda3/4.7.10
source activate pgCondaEnv

time /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_assembly_match_contigs_to_probes \
    --contigs /scratch/pg84794/UCE_run5/ext-data-genomes-transcriptomes-editnames \
    --probes /scratch/pg84794/plasm-uceprobes.fasta \
    --output /scratch/pg84794/UCE_run5/ext-data-genomes-transcriptomes-editnames-uce-search-65 \
    --keep-duplicates /scratch/pg84794/UCE_run5/ext-data-genomes-transcriptomes-editnames-duplicates65.txt \
    --min-coverage 65 \
    --min-identity 65
