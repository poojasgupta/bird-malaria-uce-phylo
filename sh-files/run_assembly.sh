#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=abyss-assembly2
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --constraint=EPYC             # node feature
#SBATCH --cpus-per-task=24             # Number of CPU cores per task
#SBATCH --mem=400G
#SBATCH --time=100:00:00
#SBATCH --output=%x_%j.out       # Standard output log
#SBATCH --error=%x_%j.err        # Standard error log

#SBATCH --mail-user=pg84794@uga.edu
#SBATCH --mail-type=NONE

cd $SLURM_SUBMIT_DIR

#module load Miniconda3/4.7.10
#source activate pgCondaEnv
#module load ABySS/2.1.5-foss-2019b
module load phyluce/1.6.8

time phyluce_assembly_assemblo_abyss \
     --config /scratch/pg84794/UCE_run5/conf-files/abyss-assembly.conf \
     --output /scratch/pg84794/UCE_run5/abyss-assembly-k55-3 \
     --kmer 55 \
     --cores 20 \
     --clean
