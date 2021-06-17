#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=match-counts-subset2-extract-align-fasta
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

#get match counts	
time  /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_assembly_get_match_counts \
    --locus-db /scratch/pg84794/UCE_run5/final-assemblies-65-abybam-abyraw-abydocap-edinames-uce-search/probe.matches.sqlite \
    --taxon-list-config /scratch/pg84794/UCE_run5/conf-files/taxon-set.conf \
    --taxon-group 'subset2' \
    --output /scratch/pg84794/UCE_run5/taxon-sets/subset2/subset2-taxa-incomp.conf \
    --extend-locus-db /scratch/pg84794/UCE_run5/ext-data-genomes-transcriptomes-editnames-uce-search-65/probe.matches.sqlite \
    --incomplete-matrix

#extract fasta
time /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_assembly_get_fastas_from_match_counts \
    --contigs /scratch/pg84794/UCE_run5/final-assemblies-65-abybam-abyraw-abydocap-edinames \
    --locus-db /scratch/pg84794/UCE_run5/final-assemblies-65-abybam-abyraw-abydocap-edinames-uce-search/probe.matches.sqlite \
    --match-count-output /scratch/pg84794/UCE_run5/taxon-sets/subset2/subset2-taxa-incomp.conf \
    --output /scratch/pg84794/UCE_run5/taxon-sets/subset2/subset2-incomplete.fasta \
    --incomplete-matrix /scratch/pg84794/UCE_run5/taxon-sets/subset2/subset2-incomplete.incomplete \
    --extend-locus-db /scratch/pg84794/UCE_run5/ext-data-genomes-transcriptomes-editnames-uce-search-65/probe.matches.sqlite \
    --extend-locus-contigs /scratch/pg84794/UCE_run5/ext-data-genomes-transcriptomes-editnames

#align fasta
time /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_align_seqcap_align \
    --fasta /scratch/pg84794/UCE_run5/taxon-sets/subset2/subset2-incomplete.fasta \
    --output /scratch/pg84794/UCE_run5/taxon-sets/subset2/subset2-mafft-nex-internal-trimmed \
    --taxa 31 \
    --aligner mafft \
    --cores 16 \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --max-divergence 0.4
