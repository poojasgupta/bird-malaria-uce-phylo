#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=bwa-host-mapping-dilcap-unen
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --constraint=EPYC             # node feature
#SBATCH --cpus-per-task=16             # Number of CPU cores per task
#SBATCH --mem=400G
#SBATCH --time=100:00:00
#SBATCH --output=%x_%j.out       # Standard output log
#SBATCH --error=%x_%j.err        # Standard error log

#SBATCH --mail-user=pg84794@uga.edu
#SBATCH --mail-type=NONE

cd $SLURM_SUBMIT_DIR

module load Miniconda3/4.7.10
source activate pgCondaEnv

time /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_snp_bwa_multiple_align \
    --config /scratch/pg84794/UCE_run5/conf-files/bwa-host-mapping-dilc-unen.conf \
    --output /scratch/pg84794/UCE_run5/bwa-host-mapping-dilcap-unen \
    --subfolder split-adapter-quality-trimmed \
    --mem \
    --cores 12
