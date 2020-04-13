Annotation Enhancer
====================
# Description
HyPAR (Hypothetical Protein Annotation Reviser) program was written to identify wrongly annotated hypothetical proteins.
The idea was based on the paper by Lozano et al. (Reference) were the authors found that approximately 70 genes that were annotated as functional proteins in one stain, were annotated as one or two smaller ORFs coding for hypothetical proteins on another strain. 
The program aims at finding frameshift variations on genes of the query genome compared to the genomic sequences of at least 5 related bacterial strains (used as reference). Then it uses blast+ to search for those sequences on the query genome and asses if they hit with one or more ORFs, and also compares the size and annotations of the gene.

## Interpretation of the resutls
* Frameshift variants are detected between the reference genome and a series of genomes belonging to the same species. A multifasta file is generated which contains the CDS, of all the analyzed strains, which presented a frameshift variant on the reference strain.

(Figure)[link]

* Local Blast searches are performed using the previously generated fasta file, against the CDS from the reference genome (.fna) and against the genomic sequence.

** There are several case scenrarios:

(Figure2)[link]


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
