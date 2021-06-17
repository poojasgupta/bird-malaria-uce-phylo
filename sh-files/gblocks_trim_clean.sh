#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=gblocks-trim-clean-all
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

#module load phyluce/1.6.8	
module load Miniconda3/4.7.10
source activate pgCondaEnv

time /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments /scratch/pg84794/UCE_run5/taxon-sets/all/all-mafft-nex-internal-trimmed \
    --b2 0.5 \
    --output /scratch/pg84794/UCE_run5/taxon-sets/all/all-mafft-nex-gblocks \
    --cores 16

time /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_align_get_align_summary_data \
    --alignments /scratch/pg84794/UCE_run5/taxon-sets/all/all-mafft-nex-gblocks \
    --show-taxon-counts \
    --cores 16

time /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments /scratch/pg84794/UCE_run5/taxon-sets/all/all-mafft-nex-gblocks \
    --output /scratch/pg84794/UCE_run5/taxon-sets/all/all-mafft-nex-gblocks-clean \
    --cores 16