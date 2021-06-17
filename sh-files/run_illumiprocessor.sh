#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=illumiprocessor-unen
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --constraint=EPYC             # node feature
#SBATCH --cpus-per-task=24             # Number of CPU cores per task
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out       # Standard output log
#SBATCH --error=%x_%j.err        # Standard error log

#SBATCH --mail-user=pg84794@uga.edu
#SBATCH --mail-type=NONE

cd $SLURM_SUBMIT_DIR

module load phyluce/1.6.8
module load illumiprocessor/2.0.9
#module load Miniconda3/4.7.10
#source activate pgCondaEnv

time illumiprocessor \
    --input /scratch/pg84794/UCE_run5/UCE_data/Mar2020_gz \
    --output /scratch/pg84794/UCE_run5/clean-fastq-unen \
    --config /scratch/pg84794/UCE_run5/conf-files/illumiprocessor_unen.conf \
    --cores 16 \
    --r1-pattern .1.fq \
    --r2-pattern .2.fq
