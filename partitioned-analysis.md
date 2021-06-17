# Partitioned Analysis

This step can be performed after sequence alignment using MAFFT and internal trimming using GBlocks. Make sure the multiple sequence alignments are in fasta or phylip format.
If not, convert nexus alignments to fasta format using the following script
```bash
phyluce_align_convert_one_align_to_another \
    --alignments /scratch/pg84794/UCE_run5/taxon-sets/subset1/subset1-mafft-nex-gblocks-clean-75p \
    --output /scratch/pg84794/UCE_run5/taxon-sets/subset1/subset1-mafft-fas-gblocks-clean-75p \
    --input-format nexus \
    --output-format fasta
```

Next, we will concatenate multiple sequence alignment using [**FASconCAT-G**](https://github.com/PatrickKueck/FASconCAT-G). This program accepts fasta or phylip format. First copy 'FASconCAT-G_v1.04.pl' to the directory of fasta sequence alignments, then run the following script
```bash
perl FASconCAT-G_v1.04.pl -s -l -j -p -p -a
```
This will output a concatenated Supermatrix in phylip format (relaxed, use just -p for strict phylip) along with a RAxML partition file (-l option) and summary of parsimony informative sites of supermatrix sequences (-j option)
		
```bash	
perl FASconCAT-G_v1.04.pl -s -j -n -n  -l -a
```
This will output nexus format (as for MrBayes) of concatenated sequence alignment. I slightly edited the partition file created for MrBayes *subset1-75p-supermatrix-partition.txt* (see example below) and nexus output of concatenated sequence alignment *subset1-75p-supermatrix.nex* for downstream processing.
```
[loci]
charset uce-1990.fasta = 26221-26602;
charset uce-295.fasta = 36328-36626;
charset uce-2105.fasta = 27489-28170;
charset uce-2213.fasta = 29777-30465;
charset uce-1800.fasta = 20876-21143;
charset uce-1588.fasta = 14366-14846;
charset uce-367.fasta = 38217-39061;
charset uce-1753.fasta = 19011-19609;
charset uce-1059.fasta = 4087-4909;
charset uce-2286.fasta = 31490-32146;
charset uce-1919.fasta = 24680-25054;
charset uce-184.fasta = 21604-21935;
charset uce-492.fasta = 42589-43474;
charset uce-805.fasta = 49613-50130;
charset uce-363.fasta = 37606-38216;
charset uce-2099.fasta = 26966-27488;
charset uce-1843.fasta = 21936-23014;
charset uce-1607.fasta = 14847-15543;
```

Now, we are ready to perform a partitioned analysis. We will use two approaches for this. First, we conducted a codon based partitioned analysis in [**PartitionFinder2**](https://github.com/brettc/partitionfinder) to find the best-fit partitioning scheme and input it to RAxML. Second, we performed an automatic selection of partitions through the SWSC-EN method, as implemented in [**PFinderUCE-SWSC-EN**](https://github.com/Tagliacollo/PFinderUCE-SWSC-EN). This program attempts to split each UCE loci into 3 parts - a core (highly conserved) and two flanking regions (generally more variable). The output from this program is used as an input to [**PartitionFinder2**](https://github.com/brettc/partitionfinder), that optimizes the partitioning scheme and can be further used for subsequent phylogenetic analysis in RAxML.

Here, I have outlined the steps for [**PFinderUCE-SWSC-EN**](https://github.com/Tagliacollo/PFinderUCE-SWSC-EN) and [**PartitionFinder2**](https://github.com/brettc/partitionfinder)

To create an input for **PFinderUCE-SWSC-EN**, we will merge the concatenated nexus alignment with the partition file. We also need to make sure that there are no -hyphens in the loci name in the partitioned section. The program expects underscores here.
```
#PFinderUCE-SWSC-EN
#NEXUS

begin data;
dimensions ntax= 38 nchar= 54622;
format datatype=dna interleave missing=-;
matrix
P_cyn_s0  -------------------- -------------------- -------------------- -------------------- --------------------
DOT1      NNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNN
DOT2      AACACGTGGACAAATTAAGA TTACAAATTATTTCAATTGT TTTTTCTTTTCTTCTACTAC CCATTCATCTGTGTGTGTGA AAAAATCTAAAGCTTCTTTA
DOT3      AACATGTAGATAAATTAAGA TTACAAATAACTTCAATTTC TTTTTTTTCTTTTCTACTGC CCATTCATCCGTGTGCGTAA AAAAATCTAAAGCCTCTTTA
DOT4      AACATGTAGATAAATTAAGA TTACAAATAACTTCAATTTC TTTTTTTTCTTTTCTACTGC CCATTCATCCGTGTGCGTAA AAAAATCTAAAGCCTCTTTA
MSBB5     AACAAGTAGACAAATTCAGA TTACAGATGATTTCAATTTT TTTTTTTTTTCTTCTACTAC CCATTCATCTGCGTGTGTAA AAAAATCTAATGCTTCTTTG
H_DB09    NNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNN
P_homo    AACATGTAGATAAATTTAAA TTACATATAATTTCTATTTT TTTTTTTTCTCAAGTACTTC CCATTCTTCAGCATGCATAA AAAAATCAAGTGCTTCTTTA
DOT5      NNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNN
H_tat     AACACGTGGACAAATTAAGA TTACAAATTATTTCAATTGT TTTTTCTTTTCTTCTACTAC CCATTCATCTGTGTGTGTGA AAAAATCTAAAGCTTCTTTA
MSB68     -------------------- -------------------- -------------------- -------------------- --------------------
P_rel_s2  AACAAGTAGCTAAATTTAAA TGACATATAATTTCTATTTT TTTTTTTTTTCGAGTACTTC CCATTCTTCTGCATGCATAA AGAAATCAAGAGCTTCTTTA
P_kno_s0  -------------------- -------------------- -------------------- -------------------- --------------------
...

;
end;

begin sets;

[loci]
charset uce_1990 = 26221-26602;
charset uce_295 = 36328-36626;
charset uce_2105 = 27489-28170;
charset uce_2213 = 29777-30465;
charset uce_1800 = 20876-21143;
charset uce_1588 = 14366-14846;
charset uce_367 = 38217-39061;

...

end;
```

This program requires Python 3.6.x or higher and several dependencies so make sure, all have been installed correctly. More details can be found on the program github page [**PFinderUCE-SWSC-EN**](https://github.com/Tagliacollo/PFinderUCE-SWSC-EN). In my case, I created a new conda environment *SWSC-EN* for installing all the necessary programs.

```bash
source activate SWSC-EN
python SWSCEN.py /scratch/pg84794/UCE_run5/PFinderUCE-SWSC-EN-master/py_script/subset1-75p-supermatrix.nex
```
This program is equipped with a nice progress bar and took about 30 mins when run directly on command line. This will genreate an output .cfg which can then be used directly as an input for [**PartitionFinder2**](https://github.com/brettc/partitionfinder) which will result in the best-fit partitioning scheme to be used for RAxML partitioned analysis.

Run **PartitionFinder2** with the *subset1-75p-supermatrix.phy* and *partition_finder.cfg* file using the following script. The input files and the script are in the same folder (in my case they are in 'partition-finder').

```bash
module load PartitionFinder/2.1.1-foss-2019b-Python-2.7.16
module load RAxML/8.2.11-foss-2019b-pthreads-sse

PartitionFinder.py /scratch/pg84794/UCE_run5/partitioned-analysis/partition-finder --raxml -p 20
```
The program will output all the result files in a separate 'analysis' folder. Within this folder, you will find *best_scheme.txt* which contains the best-fit partitioning scheme to be used as input for different tree building programs (RAxML, MrBayes, IQtree).



