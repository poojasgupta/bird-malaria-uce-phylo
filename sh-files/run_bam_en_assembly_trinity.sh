#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=trinity-bam-en-assembly
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

module load Miniconda3/4.7.10
source activate pgCondaEnv
module load Trinity/2.10.0-foss-2019b-Python-3.7.4
#module load phyluce/1.6.8

time /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_assembly_assemblo_trinity \
     --config /scratch/pg84794/UCE_run5/conf-files/cl-bam-assembly.conf \
     --output /scratch/pg84794/UCE_run5/trinity-bam-en-assembly \
     --cores 20 \
     --clean

