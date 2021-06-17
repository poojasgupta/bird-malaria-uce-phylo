# Bird-malaria-uce

This repository contains a series of commands and scripts I used for processing sequence capture data for malaria parasites. The sequence capture data was generated using custom malaria UCE probes targeting conserved loci across all representative Plasmodium and Haemoproteus (malaria parasite) genomes. My goal was to obtain good quality genomic data for avian malaria parasites and resolve its evolutionary relationships with other mammalian malaria parasites. I have mostly followed steps from Phyluce package (Faircloth 2015) but includes some additional scripts that are not available from Phyluce. The analysis described in this repository is still in progress and some of the scripts are specific to the environment that the analysis is being run on. 

Here is the bioinformatics workflow that I used -
![UCE-workflow-3](https://user-images.githubusercontent.com/55209373/122352278-d62c2180-cf03-11eb-8661-12b96b99f6eb.jpeg)

## Raw Read Processing

1a) Count number of Reads1 in raw fastq files
If the file is named as CU_084.R1.fq.gz, use:

   For a single unzipped file:

```bash
echo $(cat CU_084.R1.fq|wc -l)/4|bc
```

   For multiple files, use a loop:

```bash
for i in *.1*.fq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done
```

1b) To merge fastq files from different sequencing runs (.1 and .2 represnt Read1 and Read2). However, many forum posts suggest to merge after mapping (basically merge bam files using samtools) 

```bash
cat PFA01_2.1.fq  PFA01_5.1.fq > PFA01_merged.1.fq
cat PFA01_2.2.fq  PFA01_5.2.fq > PFA01_merged.2.fq
```

2) Adapters: 
   i7 adapter sequence is generally found in Read1 (included Rev.complement of i7 primer sequence) and
   i5 adapter sequence is generally found in Read2.
   The adapters used for malaria-UCE libraries are Illumina TruSeq HT (double-indexed). The tag sequences should be input in 5’ to 3’ orientation.
   
 ```bash
[adapters]
i7:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT*GTGTAGATCTCGGTGGTCGCCGTATCATT
``` 

3) Clean and trim raw reads: remove poor quality reads, trim poor quality base calls, and remove adapter sequences using illumiprocessor, which uses trimmomatic.
   Filtering options used for each Read(LEADING:5 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:25)
  * Clip read when average base quality over a 4bp sliding window drops below 15
  * Clip leading and trailing bases if base quality below 5 and 15, respectively
  * Skip read if shorter than 25bp

4) Check the quality of untrimmed and trimmed reads in **FastQC**.

```bash
$ find . -name '*.fastq.gz' | xargs -I % sh -c 'fastqc -o ./fastqc_out %'
```
For fastq files from a specific directory (in my case 'split-quality-adapter-trimmed'). Create the output directory ("split-adapt-trim-fastqc_out") if it doesnt exist. This is run in the clean-fastq folder.
```
grep split-adapter-quality-trimmed | find . -name '*.fastq.gz' | xargs -I % sh -c 'fastqc -o ./adapt-trim-fastqc-out %'   
#"|grep PFA|" for a specific sample
```
After creating FastQC results, run MultiQC within that same folder "adapt-trim-fastqc_out" to summarize FastQC results of all samples and get summary excel sheets.
```
multiqc .
```
I did notice that some of the reads still had some adapter contamination towards the end of the reads so probably Trimmotatic is not doing a very good job of removing adapter contamination. A web blog post here (https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html) suggested thatTrimGalore, a convenient wrapper for cutadapt, might be better at completely removing adapter contamination in reads. Yet to give it a try.

## Mapping Reads to Reference Genome
  
5) Align clean/trimmed reads to Plasmodium and Haemoproteus reference genome sequences using bwa-mem (default settings).
* [Sh-file](sh-files/run_bwa_mem.sh)
* [Conf-file](conf-files/bwa_mem_parmapping.conf)

6) Get summary stats using **Samtools** samstats and **Qualimap**. Results were further summarized in **MultiQc**.
Samtools summary stats (within individual bam outputs):
```bash
find . -name *M.bam | xargs -I % sh -c 'plot-bamstats -p %_plot_bamstats_ref/ %_samstats.txt'
```
Generate combined summary across all .bam files
```bash
find . -name *M.bam | xargs -I % sh -c 'printf "****************%*********************\n"; samtools stats % | grep ^SN | cut -f 2-; printf "\n\n********************************************\n\n"' > bwa_combined_summary_out.txt
```
```
multiqc .
```
For Qualimap stats:
```
find . -name *M.bam | xargs -I % sh -c 'printf "\n\n****************%*********************\n"; qualimap bamqc -bam % -outformat pdf;'
```
```
multiqc .
```

7) As we are only interested in keeping the reads that mapped to the parasite genome, we can use **Samtools** to remove unmapped reads, keep the mapped reads. The sam file outputs Read IDs in the first column and mapped/unmapped status in the 4th column - usually ‘0’ for unmapped reads and a non-zero for mapped reads. The -f/-F options to the samtools refers to the presense/absence of bits in the FLAG field. So -f 4 looks for '0' and only output alignments that are unmapped (flag 0×0004 is set) and -F 4 looks for non-zero in the 4th column and only output alignments that are not unmapped (i.e. flag 0×0004 is not set), hence these would only include mapped alignments. 

This needs to be run in the 'bwa-par-mapping' dir.
We can count the total number of mapped reads with the samtools -F 4 flag. Alternativley, we can count only the unmapped reads with -f 4

### Extract mapped Reads only:
```bash
find . -name *M.bam | xargs -I % sh -c 'samtools view -F 0x04 -b % > %-mapped.bam' 
```

### Extract unmapped Reads only:
```bash
find . -name *M.bam | xargs -I % sh -c 'samtools view -f 0x04 -b % > %-unmapped.bam' 
```

To verify results, count the mapped and unmapped reads and match with the output from Qualimap,
```bash
samtools view -c CU_084_P-CL-RG-MD-M.bam-mapped.bam
samtools view -c CU_084_P-CL-RG-MD-M.bam-unmapped.bam
```

8) The bam file containing the mapped reads needs to be sorted first before converting it back to Fastq files. We will sort the mapped.bam file by name/read group.
For a single sample,
```
samtools sort -n CU_084_P-CL-RG-MD-M.bam-mapped.bam -O bam -o CU_084_P-map-sort.bam
```
For multiple samples,
```
find . -name '*mapped.bam' | xargs -I % sh -c 'samtools sort -n % -O bam -o %-sort.bam'	
```
Rename files to something shorter,
```
find . -name '*sort.bam' |
    while IFS= read i; do
        cp "$i" "${i%%-CL*}-map-sort.bam";
    done
 ```
9) Now, we can convert these mapped-sorted bam files to fastq files by multiple tools such as Picard, bedtools, samtools. Most of these programs discard the secondary alignments and singleton reads, so the total number of reads obtained are less than those found in original mapped bam alignment file. After testing all programs, Samtools works best in outputting the number of mapped reads. It can be configured to give Reads1, Reads2, singletons and other secondary alignments, thus matching closely with the total number of mapped reads as obtained with .bam file.

For a single file, 

```
samtools fastq -t -1 16N1030_LHP-map-READ1.fq -2 16N1030_LHP-map-READ2.fq -s 16N1030_LHP-map-READ-singleton.fq -0 extra-reads.fq -N 16N1030_LHP-map-sort.bam
```
-t option retains the @RG header info from bam file and can be used to distinguish which read came from which sample (run1, run2 etc).
For multiple files, (run from the bwa_par_mapping folder, it will create an output folder directly in clean-fastq folder)

```
ls | xargs -I % sh -c 'filePathPrefix=/scratch/pg84794/UCE_run2/clean-fastq/$(echo "%" | rev | cut -d"/" -f1 | rev | cut -d"-" -f1)/samtools-bam; fileNamePrefix=$(echo "%" | rev | cut -d"/" -f1 | rev | cut -d"-" -f1);printf "\n\n****************%*********************\n"; mkdir $filePathPrefix; samtools fastq -t -1 $filePathPrefix/$fileNamePrefix-map-READ1.fq -2 $filePathPrefix/$fileNamePrefix-map-READ2.fq -s $filePathPrefix/$fileNamePrefix-map-READ-singleton.fq % -N;'
```
Alternatively, if you need to output in a separate folder, use the following code. ' rev | cut -c13- | rev' means reverse the file name (.bam file) cut 13 characters from the beginning (technically last 13 characters) and then reverse it back to generate output folder names in 'clean-bam-fastq' folder. For this to work, first create the output folder (clean-bam-fastq) in the main directory. Also make sure all .bam files are in a separate 'mapped reads folder' within the bwa_par_mapping folder.

```
ls | xargs -I % sh -c 'fileNamePrefix=$(echo "%" | rev | cut -c16- | rev); filePathPrefix=/scratch/pg84794/UCE_run5/clean-docap-bam-fastq/$fileNamePrefix;printf "\n\n****************%*********************\n"; mkdir $filePathPrefix; samtools fastq -t -1 $filePathPrefix/${fileNamePrefix}_unmap-READ1.fq -2 $filePathPrefix/${fileNamePrefix}_unmap-READ2.fq -s $filePathPrefix/${fileNamePrefix}_unmap-READ-singleton.fq -0 $filePathPrefix/${fileNamePrefix}_unmap-extra-reads.fq % -N;'

```
Move all the extra reads in a separate folder from the clean-docap-bam-fastq folder otherwise it throws an error during assembly
```
find . -name *unmap-extra-reads.fq | xargs -I % sh -c 'mv % extra-reads/.'
```

## Reads Assembly

10) We can assemble reads to contigs using two apporaches- De novo assembly and reference based assembly. Following the standard phyluce pipeline and default settings, I used trinity and Abyss (kmer=31 and 55, Try KmerGenei to get optimal Kmer for the dataset) to assemble my cleaned/trimmed raw reads as well as ran assembly with the mapped reads obtained from bam alignement file. I was interested in comparing whether mapped reads assembled better compared to all the reads containing both parasite and host reads. Additionally, I wanted to compare the two assemblers as Trinity is primarliy for RNA-seq reads while Abyss is dedicated for DNA genome assembly.
   
## Summary assembly stats
```
for i in adapt-trim-abyss-assemblies-k55/contigs/*.fasta; do phyluce_assembly_get_fasta_lengths --input $i --csv; done
```

11) Run assembly QC and get summary stats using **MultiQC**

## Finding UCE loci in assembled contigs

12) After assembling the contigs using abyss or trinity, we will search for contigs that match our UCE loci of interest. To do this, we will match contigs from each taxon to probes.fasta file (plasm-uceprobes.fasta) which contains the probe set used for the enrichments. As some taxa are distantly related, I used slightly relaxed matching parameters (set min identity and min coverage as 65).
```
> source activate pgCondaEnv
> time /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_assembly_match_contigs_to_probes \
    --contigs /scratch/pg84794/UCE_run5/assemblies/abyss-bam-en-assembly-k55/contigs \
    --probes /scratch/pg84794/plasm-uceprobes.fasta \
    --output /scratch/pg84794/UCE_run5/uce-search/abyss-bam-en-assembly-k55-uce-search \
    --keep-duplicates /scratch/pg84794/UCE_run5/uce-search/abyss-bam-en-assembly-k55-duplicates.txt \
    --min-coverage 65 \
    --min-identity 65
```

A summary of UCE loci matching contigs for each taxon can be obtained using the following utility script (This script was sourced from https://github.com/MikeWLloyd/Utility-Scripts) which can then be opened as csv. This is essentially that valuable capture data which will tell us whether our target enrichment method worked or not.

```
python match_contigs_log_parse phyluce_assembly_match_contigs_to_probes.log abyss-bam-docap-uce-search-out.txt
```
## Extracting UCE loci

13) After identifying UCE loci, we will extract the fasta sequence matching UCE loci. But first, we would like to incorporate parasite data from external sources to test the utility of our UCE probes againsts exemplar genomes and broaden taxon sampling for phylogenetic inference.

### Extracting UCE loci from external data (genomes/transcriptomes)

13a) This step allows you to incoporate genomic/transcriptomic data from external sources (e.g. avilable on NCBI, VEuPathDB, MalAvi) or from previous studies in your phylogenomic analysis. This step can also be used to incorporate outgroup data for your trees.

For my analysis, I downloaded genome and transcriptome data for haemosporidian parasites from NCBI and [MalAvi](http://130.235.244.92/Malavi/Downloads/) databases. The data available on NCBI were primarily genome assemblies (*.fna* files) which need to be converted first to a simpler 2bit format before it can be used for downstream analysis. Here I have followed steps similar to those outlined in [**phyluce tutorial 3**](https://phyluce.readthedocs.io/en/latest/tutorial-three.html#tutorial-iii-harvesting-uce-loci-from-genomes)

```
$ babesia bovis 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/165/395/GCF_000165395.1_ASM16539v1/GCF_000165395.1_ASM16539v1_genomic.fna.gz

$ Theileria annulata
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/225/GCF_000003225.3_ASM322v1/GCF_000003225.3_ASM322v1_genomic.fna.gz

$ Plasmodium reichenowi
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/601/855/GCF_001601855.1_ASM160185v1/GCF_001601855.1_ASM160185v1_genomic.fna.gz

$ Plasmodium malariae
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/090/045/GCF_900090045.1_PmUG01/GCF_900090045.1_PmUG01_genomic.fna.gz

$ Plasmodium vivax
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/415/GCF_000002415.2_ASM241v2/GCF_000002415.2_ASM241v2_genomic.fna.gz

$ Plasmodium ovale
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/090/025/GCA_900090025.2_PowCR01/GCA_900090025.2_PowCR01_genomic.fna.gz

$ Plasmodium gaboni- chimps
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/602/025/GCF_001602025.1_ASM160202v1/GCF_001602025.1_ASM160202v1_genomic.fna.gz

$ Plasmodium vinckei
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/681/995/GCF_900681995.1_PVVCY_v1/GCF_900681995.1_PVVCY_v1_genomic.fna.gz

$ Plasmodium fragile
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/956/335/GCF_000956335.1_Plas_frag_nilgiri_V1/GCF_000956335.1_Plas_frag_nilgiri_V1_genomic.fna.gz
```
First, unzip all '.gz' files and convert them to 2 bit format

```
gunzip *.gz
```
We will convert these fasta files to *2bit* format using faToTwoBit utility. You can download the utility files from the [Kent Source Archive(https://hgdownload.soe.ucsc.edu/admin/exe/)

```
./faToTwoBit P_gaboni.fna P_gaboni.2bit
```
Next, we need to make sure that each 2bit genome file is within its own individual directory

$ P_gaboni >> P_gaboni.2bit
$ P_frag >> P_frag.2bit

## Find UCE loci and extract fasta

Now, we will use *find_loci.sh* script to find uce matches across the various genomes and store the search results in lastz file format as well as summary of those search results in a *sqlite* database
```
> source activate pgCondaEnv
> /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_probe_run_multiple_lastzs_sqlite \
    --db /scratch/pg84794/genomes/parasites/ext-data2.sqlite \
    --output /scratch/pg84794/genomes/parasites/ext-data-genome2-lastz \
    --scaffoldlist P_frag P_gaboni P_mala P_ovale P_reich P_vinc P_vivax Hepatocyst P_delich P_ashfordi H_colum \
    --genome-base-path ./ \
    --probefile /scratch/pg84794/plasm-uceprobes.fasta \
    --cores 12
```

After identifying UCE loci across genomnes, we will extract FASTA sequence matching UCE loci from genome sequences along with some sequence flanking each UCE locus using the following config file and script
```
$ config file
[scaffolds]
P_frag:/scratch/pg84794/genomes/parasites/P_frag/P_frag.2bit
Hepatocyst:/scratch/pg84794/genomes/parasites/Hepatocyst/Hepatocyst.2bit
P_gaboni:/scratch/pg84794/genomes/parasites/P_gaboni/P_gaboni.2bit
P_mala:/scratch/pg84794/genomes/parasites/P_mala/P_mala.2bit
P_ovale:/scratch/pg84794/genomes/parasites/P_ovale/P_ovale.2bit
P_reich:/scratch/pg84794/genomes/parasites/P_reich/P_reich.2bit
P_vinc:/scratch/pg84794/genomes/parasites/P_vinc/P_vinc.2bit
P_vivax:/scratch/pg84794/genomes/parasites/P_vivax/P_vivax.2bit
P_delich:/scratch/pg84794/genomes/parasites/P_delich/P_delich.2bit
H_colum:/scratch/pg84794/genomes/parasites/H_colum/H_colum.2bit
```

```
$ script
> phyluce_probe_slice_sequence_from_genomes \
    --lastz /scratch/pg84794/genomes/parasites/ext-data-genome2-lastz \
    --conf /scratch/pg84794/genomes/parasites/scaffold2.conf \
    --flank 500 \
    --name-pattern "plasm-uceprobes.fasta_v_{}.lastz.clean" \
    --output /scratch/pg84794/genomes/parasites/ext-genome2-fasta
```
The extracted fasta sequence for each taxa in 'ext-genome2-fasta' directory can now be used for further downstream analyses and in my case, I combined it with other transcriptome data (sourced from MalAvi) to generate the external parasite database

I ran step 12 to obtain UCE loci matches against these reference fasta sequences (obtained in step 13a) and transcriptome dataset for avian haemosporidian parasites. This will output *lastz* files for each taxon and probes.matches.sqlite database
```
> time /home/pg84794/.conda/envs/pgCondaEnv/bin/phyluce_assembly_match_contigs_to_probes \
    --contigs /scratch/pg84794/UCE_run5/ext-data-genomes-transcriptomes \
    --probes /scratch/pg84794/plasm-uceprobes.fasta \
    --output /scratch/pg84794/UCE_run5/ext-data-genomes-transcriptomes-uce-search \
    --keep-duplicates /scratch/pg84794/UCE_run5/ext-data-genomes-transcriptomes-duplicates.txt \
    --min-coverage 65 \
    --min-identity 65
```

14) Now, in order to extract the fasta sequence matching UCE loci from our study taxa as well as external dataset, we will first create a data matrix configuration file specifying the taxa and loci to be included for subsequent analysis. We can get a good idea about which taxa to include based on the results of 'Finding UCE loci' step. I chose to exclude taxa that worked poorly and aimed to generate an 'incomplete matrix' with 'dataset1' containing taxa withe than 50 loci recovery. I also included other datasets in my data matrix configuration file (datasets containing more than 100loci recovery, datasets excluding mixed infections, etc).

## Aligned UCE loci across samples- summary stats

15) Run from the folder which contains mafft-nexus-internal-trimmed folder. 
```
phyluce_align_get_align_summary_data  \
--alignments subset2-mafft-nexus-internal-trimmed \
--output-stats subset2-alignment-stats.csv \
--show-taxon-counts
```
