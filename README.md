# FASTKs
Pipeline for generating Ks plots of transcriptome data
=============

<b>Author History</b>:
Original by Kerr Wall: Sept 21 2005 					  
Modified by Liying Cui, Severn Everett 			   
Modified by Michael R. McKain: Jun 29 2011               			   
	Removed ESTscan: May 23 2012					   
	Rewrite: Feb 16 2013				

version 1.1, October 06, 2015<br>
<b>Contact</b>: https://github.com/mrmckain

<h3>Description</h3>
FASTKs estimates pairwise Ks, Ka, and Ka/Ks values for transcriptome datasets. FASTKs can be used to estimate these statistics for paralogs (a single dataset) or orthologs (two datasets). 

<h4>Requirements</h4>

All required software and the FASTKs directory are expected to be in the user's path.

BLAST+ and <a href="http://abacus.gene.ucl.ac.uk/software/paml.html"><b>PAML</b></a>.

The program <a href="https://www.ncbi.nlm.nih.gov/pubmed/16845082"><b>PAL2NALv14</b></a> is packaged with FASTKs.

The Perl module Parallel::ForkManager.

The included R script for normal mixture modeling and producing a Ks frequency plot requires <a hre="http://www.stat.washington.edu/mclust/"><b>mclust</b></a>.


<h4>Input</h4>

<b>Preparation of Input</b>:<br>
Transcript CDS and accompaning peptide sequences should be in two files: XXX.cdna and XXX.pep. XXX can be any name, but it must be congruent between the two files. 

A blast output file must be provided. The query for the blast must be the query sequence (see below) used for the pipeline.  The database for the blast must be the subject sequence (see below). The output should be in tabular format. Here is an example command line blast run using blastn from BLAST+:

	makeblastdb -in YYY.cdna -dbtype nucl
	blastn -num_threads 2 -query XXX.cdna -db YYY.cdna -outfmt 6 -evalue 1e-40 > XXX-YYY.blastn

<h4>Output</h4>

FASTKs will create a directory named XXX_YYY_ks.  If this directory already exists, FASTKs will rename the existing directory XXX_YYY_ks.1. If XXX_YYY_ks.1 exists, FASTKs will rename XXX_YYY_ks XXX_YYY_ks.2 and any subsequent directories will be renamed in increasing numerical order.  The latest iteration of the run will alway be name XXX_YYY_ks.  

The two output files that will be of interest to the user are XXX_YYY.paralogs.txt and XXX_YYY.kaks.txt.

		XXX_YYY.paralogs.txt 	Contains the query and subject sequences in a homolog pair and their respective overlapping positions.
		XXX_YYY.kaks.txt 	Contains columns of Pair number, Ka, Ks, Ka/Ks, query sequence, and subject sequence. This file is used the in the run_mclust.r script to produce the Ks frequency plot.


<h4>Usage</h4>

perl FAST_kspipe.pl --query CDS_file --subject CDS_file --blast blast_output [options]

<b>Options:</b>
         
	-query          First transcriptome dataset for Ks plot estimation. 
	-subject        Second transcriptome dataset for Ks plot estimation.  If the user wants to estimate Ks for paralogs, the query file should also be used here. 
	-blast_output   Blastn output file in tabular format. Files should be make using blastn and the -outfmt 6 option. The query file (above) should be used to query a database made from the subject file (above).
	-pairs          Total number of unique pairs to be included in sets that are passed to ForkManager. The more processors you use the lower this number can be, and the faster FASTKs will finish. [Default = 500]
	-perid          Minimum percent identify acceptable to pass a putative pair. [Default = 40]
	-alen           Minimum alignment length with a transcipt pair. [Default = 100]
	-processes      Total number of processors to be used for ForkManager. [Default = 2]

<h4>Other Scripts</h4>

<b>Batch Paralog Ks Frequency Plot Estimation</b>

<h5>Usage</h5>
./run_kspipeline.sh

	Script must be run from a directory that has CDS and peptide files for transcriptomes. Files should be in the format AAA.cdna and AAA.pep. Any number of these files can be in the directory.  This script will create subdirectories for each pair, run blastn (expects BLAST+ to be in path), and then run FASTKs.

<b>R Function for Estimating Normal Mixtures and Plotting with Ks Frequency Plot</b>

<h5>Usage</h5>

Load script run_mclust.r in R.  If user wants to run command line, additions need to be made. 

runmclust(ksfile,comp,species)
	
	-ksfile          Path to XXX_YYY.kaks.txt output from FASTKs.
	-comp			 Number of components for mixture modeling. Recommendation: [4]
	-species		 Taxon name to be used as main title of plot.

<h4>Citing</h4>

FASTKs can be cited through this repository or <a href="http://gbe.oxfordjournals.org/content/8/4/1150.long"><b>McKain et al. 2016. A Phylogenomic Assessment of Ancient Polyploidy and Genome Evolution across the Poales. Genome Biology and Evolution 8(4):1150-1164.</b></a>

In addition, the following should be cited:

<a href="http://nar.oxfordjournals.org/content/34/suppl_2/W609.long"><b>Suyama M, Torrents D, Bork P. 2006. PAL2NAL: robust conversion of
protein sequence alignments into the corresponding codon alignments.
Nucleic Acids Res. 34:W609–W612.<b></a>

<a href="http://mbe.oxfordjournals.org/content/24/8/1586.abstract"><b>Yang Z. 2007. PAML 4: phylogenetic analysis by maximum likelihood. Mol Biol Evol. 24:1586–1591.</b></a>

If using the R script for normal mixture modeling, cite:

<a href="https://www.stat.washington.edu/raftery/Research/PDF/fraley2002.pdf"><b>Fraley C, Raftery AE. 2002. Model-based clustering, discriminant analysis, and density estimation. J Am Stat Assoc. 97:611–631.</b></a>

<a href="https://www.stat.washington.edu/research/reports/2012/tr597.pdf"><b>Fraley C, Raftery AE, Murphy TB, Scrucca L. 2012. mclust Version 4 for R: normal mixture modeling for model-based clustering, classification, and density estimation.</b></a>
	