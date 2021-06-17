### Find all files with '.bam' and run 'samtools flagstat' on each of those files (Denoted by %) and output them to '.txt'

```
find . -name *.bam | xargs -I % sh -c 'samtools flagstat % > %_flagstat.txt'

```

### Renaming files
basename is a function in UNIX that is helpful for removing a uniform part of a name from a list of files.
For example, we can use basename to remove the .fastq.gz extension from the files.

basename CU_084_P.fastq.gz .fastq.gz
We see that this returns just the SRR accession, and no longer has the .fastq file extension on it.

CU_084

### Search and replace in vim:

```
:%s/<term-to-replace>/<what-to-replace-with>/g
```

/g means replace globally without confirmation
/gc it'll ask for confirmation for each replacement

### NGS Analysis

### Command for displaying file names in a directory. Could be used for creating desired config files

```
ls | xargs -I % sh -c 'printf "\n\n****************%*********************\n";'
```

### Command used for generating assembly.conf file from folders containing clean 'bam to fastq' reads
```
ls | xargs -I % sh -c 'printf "%:/scratch/pg84794/UCE_run3/clean-bam-fastq/%/\n"'
```
### Renaming files. Remove '.' and replace with '_' (For file names beginning with H; e.g., H.DH18_5_mgd.contigs.fasta)
```
ls H.* | xargs -I % sh -c 'fileNamePrefix=$(echo "%" | cut -c3-); cp % H_$fileNamePrefix;'
```

### Check if your output bam is sorted or not
```
samtools view -H IDLmerged.bam | grep "SO:"
```

### Extract @RG header info to a file
```
samtools view -H $P.IDL18-1609_1_r1_M.bam-mapped.bam | grep "^@RG" >RG_$P.IDL18-1609_1_r1
```

### Extracting all .bam files from individual folder to a separate folder (e.g. bam-out)
```
find . -name *M.bam | xargs -I % sh -c 'cp % bam-out/.'
```

### Running bamqc on this set of .bam files when extracted to a separate folder (e.g. bam-out) and getting muultiqc summary
```
find . -name '*M.bam' | xargs -I % sh -c 'printf "\n\n****************%*********************\n"; qualimap bamqc -bam % -outformat pdf;'
multiqc .
```
### Convert one alignment to another
```
phyluce_align_convert_one_align_to_another \
    --alignments /scratch/pg84794/UCE_run5/taxon-sets/subset1/subset1-mafft-nex-gblocks-clean-75p \
    --output /scratch/pg84794/UCE_run5/taxon-sets/subset1/subset1-mafft-fas-gblocks-clean-75p \
    --input-format nexus \
    --output-format fasta
```
