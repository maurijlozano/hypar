Hypothetical Protein Annotation Reviser
====================
# Description
HyPAR (Hypothetical Protein Annotation Reviser) program was written to identify wrongly annotated hypothetical proteins.
The idea was based on the paper by Lozano et al. (Not published yet) were the authors found that approximately 70 genes that were annotated as functional proteins in one strain, were annotated as one or two smaller ORFs coding for hypothetical proteins on another strain. 
The program aims at finding frameshift variations on genes of the query genome compared to the genomic sequences of at least 5 related bacterial strains (used as reference). Then it uses blast+ to search for those sequences on the query genome and asses if they hit with one or more ORFs, and also compares the size and annotations of the gene.

## Interpretation of the resutls
* Frameshift variants are detected between the reference genome and a series of genomes belonging to the same species. A multifasta file is generated which contains the CDS, of all the analyzed strains, which presented a frameshift variant on the reference strain.
* Local Blast searches are performed using the previously generated fasta file, against the CDS from the reference genome (.fna) and against the genomic sequence.
** There are several case scenarios:

![alt text](https://github.com/maurijlozano/hypar/blob/master/fig.png)

# Installation

## Required Dependencies
In order to run, HyPAR.py requires the following dependencies:  
* [Snippy](https://github.com/tseemann/snippy)
* Biopython package
* Blast+

# Running HyPAR#
usage: 
```
         HyPAR.py [-h] -G QUERYGENUS -S QUERYSPECIES [-a ACC [ACC ...]] [-f FILES [FILES ...]] -s
         QUERYSTRAIN -e YOUREMAIL [-o OUTPUT] [-n NSTRAINS]
```

The following arguments are required: -G/--Genus, -s/--Strain, -e/--email  
Folder names containing blank spaces must be surrounded by quotation marks.  
-n  argument is optional. By default the program will attempt to download the genome of 5 strains to compare with the query sequence.  
Additionally, Accession numbers or genbank and fasta files can be supplied with the -a or -f parameters. If -f parameter is active, accession numbers and web search for the query genome will be omitted.
When running with -f parameter, the program requires both genebank file and a fasta file, each with all the sequences for the bacteria (i.e. a multigenebank and multifasta files with all replicons). The files must have the same basename: e.g. genome.gb and genome.fasta

e.g.  
Run with web search only  
`./HyPAR.py -G genus -S species -s strain -e email -o ResultsFolder`
Run with query accession numbers
`./HyPAR.py -G genus -S species -s strain -e email -o ResultsFolder -a NC_003037.1 NC_003047.1 NC_003078.1`
Run with query gb and fasta files
`./HyPAR.py -G genus -S species -s strain -e email -o ResultsFolder -f name.gb name.fasta`

# Paper
Unfortunately no paper would be available since in the last few year the NCBI began updating the annotation of all the genomes, making this program rather useless. At the begining of project the program detected, on several bacterial species, more than a hundred genes which annotation could be corrected. During the last month, the program only was capalbe of detecting a few to none genes.

## INTRODUCTION
During the last decades, with the development of second and third generation sequencing technologies, raw genomic sequence data has increased exponentially reaching in February 2020 almost 400 billion bases ([NCNI Statistics](https://www.ncbi.nlm.nih.gov/genbank/statistics/)), 33 million RefSeq nucleotide accessions, and a 100 thousand complete genomes of which 60 thousand correspond to bacteria ([NCBI Refseq statistics](https://www.ncbi.nlm.nih.gov/refseq/statistics/)).
All these prokaryotic genomes were annotated, in general, using automated platforms which perform a set of steps involving prediction of RNA genes, repetitive regions (CRISPRs) and mobile genetic elements, finding protein coding sequences (CDS) and pseudogenes, and finally, assigning a predicted function to all the CDS ([Tatusova et al. 2016](https://academic.oup.com/nar/article/44/14/6614/2468204), [Galens et al. 2011](https://environmentalmicrobiome.biomedcentral.com/articles/10.4056/sigs.1223234), and others).   
Although Bacterial and archaeal genomes usually lack introns, which make the identification of genes boundaries easier, gene calling procedures are not a fail proof ([Poptsova et al. 2010](https://www.microbiologyresearch.org/content/journal/micro/10.1099/mic.0.033811-0)). Errors occur for several reasons such as the use of alternative start codons, different possible translation initiation sites (TIS)([Zhu et al. 2004](https://academic.oup.com/bioinformatics/article/20/18/3308/201914), [Yada et al. 2001](https://academic.oup.com/dnaresearch/article/8/3/97/376303)), alternative functions of stop codons -which can code for selenocystein and pyrrolysine- ([Theil Have et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3639795/)), the presence of frameshift mutations (real or arising from sequencing errors) and the limitations of HMM-based gene prediction algorithms. 

Common problems arising from this errors concern the accuracy of translation initiation sites (TIS) and the erroneous identification of short open reading frames (ORFs). In the last case, a current practice that has been used to avoid overconfidence in the annotation is the use of ‘hypothetical protein’ or ‘hypothetical conserved protein’.

Incorrectly annotated short ORFs as hypothetical conserved proteins, might be even included in HMM profiles, propagating the error to newly sequenced genomes by the automated genome annotation pipelines. Additionally, short ORFs might be annotated based only on the presence of a conserved protein domain.

Even though during the last decade the annotation pipelines have evolved to correct some of this errors - including pseudo gene annotation- ([Tatusova et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5001611/), [Tanizawa et al. 2018](https://academic.oup.com/bioinformatics/article/34/6/1037/4587587)), and the process of reannotation of genomes has begun ([Haft et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5753331/)), during our work with genomes of model rhizobia (not published) we found short proteins, mostly annotated as hypothetical, that presented partial similarity to larger genes with a more complete functional annotation on closely related strains, and for which a single nucleotide polymorphisms resulting in a frameshift had been detected. Such ORFs are correctly annotated as pseudogenes in some cases, more frequently on recently annotated genomes, which is not the case of most model bacterial strains. Correction of annotation errors like the described before, could help understand the differential behaviour of related bacterial strains.
 
Based on these findings, we developed a program, HyPAR (Hypothetical Protein Annotation Reviser) to determine whether a short ORF corresponds to an hypothetical protein or should instead be annotated as a larger pseudogene.

## Material and methods
HyPAR is a command line program implemented in phyton which uses Biopython, Pandas and numpy modules, in accordance with Entrez Programming Utilities (E-utilities), blast+ suite, and Snippy program to: download the genome of query and at least 5 reference bacterial strains (unless provided by the user, see program parameters), run Snippy in contig mode to search for single nucleotide insertions and deletions resulting in frameshift mutations in the query genome (makes all the pairwise analysis between the query and each reference genome), and make blast searches to find the gene coordinates in the query genome. As a result, a table is generated containing the identified hypothetical protein on the query genome, and the suggested correction based on all the reference genomes used.
## HyPAR algorithm
A detailed scheme of HyPAR algorithm is displayed on figure 1. As the first step HyPAR retrieves and copies the query and reference genomes in the designated folder (-o parameter). There are three different modes to input the required genomes, the first is simply to provide to HyPAR the genus (-G genus), specie (-S specie) and strain (-s strain) of the query. In this mode, HyPAR will try to download from the NCBI Assembly database the genome in genbank format for both the query and reference genomes (5 by default, -n number). As an alternative, NCBI accessions can be provided (-a accession). The third mode is to directly input the path to the already downloaded query and reference genomes (-f files).
Once the query and reference genomes are downloaded, CDS sequences and positions are extracted from the Query genbank file, and a Snippy analysis comparing the query to each reference genome in a pairwise fashion is performed. As a result a snps.csv is generated for each snippy run. The next step looks for reference genome CDSs for which a variant with the description of frameshift or stop gained has been discovered. For those cases the nucleotide sequence, locus_tag and product features are extracted from the reference genbank file, and a fasta file and a locus_tag/product table are generated. Here only the polymorphisms related to non hypothetical proteins are retained. The fasta file is then used on blastn searches against the query fasta nucleotide file containing all the CDSs (fna file used to make blast database) and the query genomic sequence (query.fasta used to make blast database). The first blast search is to gather information about the annotated genes on the query and reference genomes. The second blast search is used to search for the gene coordinates on the query genome, and to discard meaningless results, for example low query cover regions, where the sequence from the reference is not sufficiently similar to the query sequence. Several case scenarios may occur: first, the gene on the reference genome matches with a low query coverage (see Figure) to one or two genes on the query; second, very low query coverage, in the case of a missing gene and an overlapping gene; last, cases with high query coverage, which might involve genes with truncated/alternative starts codons, or with truncated ends.




