![Python 3.7](https://img.shields.io/badge/Python-3.7-blue.svg) ![Ubuntu 18.04.2](https://img.shields.io/badge/Ubuntu-18.04.2-green.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Systematic analysis of the underlying genomic architecture for transcriptional-translational coupling in prokaryotes


**Systematic analysis of the underlying genomic architecture for transcriptional-translational coupling in prokaryotes**  
Richa Bharti, Daniel Siebert, Bastian Blombach, Dominik G. Grimm

 <p style='text-align: justify;'> The aim of this study is to provide a comprehensive workflow to systematically investigate bacterial genomes for the abundance of transcriptional and translational associated genes clustered in distinct operons.</p>

 ## Pipeline Summary
We have created a comparative genomics pipeline for screening genomic distributions of probable conserved operonic motifs in bacteria (Figure 1). More details can be found in the accompanying manuscript (currently under preparation). The workflow is based an a series of different steps, which are based on custom Python, R and Bash scripts. The full pipeline can be run with a single Bash command: `run.sh` (more details can be found below).

Figure 1 gives a general overview about the different steps of the pipeline. More details can be found in the accompanying manuscript.


<p align="center">
  <img src="https://github.com/grimmlab/transcriptional-translational-coupling/blob/master/Figure1.png">
</p>

# Prerequisites
Basic prerequisites that need to be satisfied:

### OS
Any Linux based distro should work. We tested the scripts using:

Distributor ID: Ubuntu <br/>
Description:    Ubuntu 18.04.2 LTS <br/>
Release:        18.04 <br/>
Codename:       bionic <br/>

`lsb_release -a` on a Ubuntu based system.

###  Dependencies, Packages & Installation
<p style='text-align: justify;'> In order to reproduce the results please make sure that you have all the libaries on your local machine.
Without these the workflow will fail to run. </p>

1. Clone this project:
```
git clone https://github.com/grimmlab/transcriptional-translational-coupling.git
```
   
2. If not already installed, please install R (>4.0.3) and the taxize library on your local machine:
```
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo apt update
gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | sudo apt-key add -

sudo apt install r-base r-base-core r-recommended r-base-dev

sudo Rscript -e 'install.packages("taxize",dependencies = TRUE)'
```
   
3. Install all Python dependencies (we recommend to first setup a virtual Python environment): 
```
pip3 install -r requirements.txt
```
    

# Data
To identify and investigate operons containing both transcriptional and translational genes, 2,071 bacterial genomes were downloaded from the DOOR2 database and corresponding annotation files were retrieved from GenBank using the available REST API.  

We created a data dump, including all necessary data from the DOOR2 database and NCBI Genbank to reproduce the results from the paper. Alternatively all data can be also fetched from the DOOR2 database and from NCBI Genbank (time consuming).

## Download data dump
To download the data dump just run the following commands in your command line.  

**This is only needed if you do not run the full pipeline. To run the full pipeline have a look at the next section**

1. First clone this repository:
```
git clone https://github.com/grimmlab/transcriptional-translational-coupling.git
```
   
2. Move to the github and data folder:
```
cd transcriptional-translational-coupling/data
```
   
3. Merge splitted data zip files:
```
zip -FF data_dump.zip --out data.zip
```
   
4. Unzip merged zip file to extract the data
```
unzip data.zip
```

# How to run the pipeline

The full pipeline can be run by executing a single bash script (running the full pipeline can **take several hours**):

```
sh run.sh
```

This will run the full pipeline, as illustrated in Figure 1B.  


### Summary of the different steps of the pipeline
In the following we will give some detailes about the individual steps of the pipeline within the `run.sh` bash script. More details about the general pipeline can be found in the accompanying manuscript (currently under preparation):

1. The environment and all PATH variables are setup by the script.  

2. `run_uncompress_door2_and_ncbi_data`: The data dumps from the GitHub repo are merged and unzipped into `data/data`  

3. `run_classification_code`: This is the most time intensive step of the whole pipeline. Here the data from GenBank and the DOOR database are analysed. First, GenBank annotation files are extracted for each genome and a table of operons is created based on the relative proportions of operons which fall into one of the following categories:  
*a) Genes in operon associated with only transcription  
b) Genes in operon associated with only translation  
c) both: Genes in operon associated with both transcription and translation   
d) Genes in operon associated with neither transcription and nor translation.*  
Second, count data is generated for each operon table by comparing them with a list of bacterial transcriptional and translational genes. The resulting gene list consists of gene names and their reported synonyms for each individual entry. A simultaneous keyword (gene name) and synonym-based (gene-synonym) search module is utilized to create a count table containing a catalogue of each of the categories.  
6. `run_concatenate_all_classified_files`: The output files from the previous task are mergerd into a single table, including information such as locus tag, function, gene name and COG id. In addition the table is filtered for operons with only genes associated with both transcription and translation (results can be found in `analyses/genome_list_containing_both.txt` and `analyses/Final_combined_files.txt`)
7. `run_modify_concatenate_all_classified_files`: Due to various inconsistencies how genes are named and identified as well as inconsitent or missing annotations an additional script is used create a unified and cleaned up table (output can be found in `analyses/Final_combined_files_corrected.txt`) 
8. `run_occurrence_based_ranking`: In this step a occurance based ranking is performed in which genes, functions and COG ids are grouped and clustered based on their occurrence in the genomes. The top 18 overlapping occurrence genes were extracted and used to perform a gene enrichment and clustering-based analyses. Results of this analysis, including plots, can be found in `analyses/occurrence`.
9. `run_cooccurrence_based_gene_ranking`: All cooccurrencing gene motifs are grouped and counted. Results of this analysis, can be found in `analyses/co-occurrence/gene_cooccurrence_based_ranking_output.txt`.
10. `run_cooccurrence_based_functional_ranking`: All cooccurrencing functional motifs are grouped and counted. Results of this analysis, can be found in `analyses/co-occurrence/functional_cooccurrence_based_ranking_output.txt`.
11. `run_gene_motif_search`: The STRING v10 database is used to perform a network based analysis for clustering motifs based on gene fusion (genes reportedly existing as hybrids without any intergenic sequence(s)), gene neighborhood (genes within close proximity) and gene co-occurrence (genes existing together on same genomic loci with intergenic sequences and/or other genes). Next, the frequency of the resulting operonic gene motifs across all extracted genomes are computed. Results of this analysis, can be found in `analyses/analyses/motif_counts_for_genes.txt` and in `analyses/gene_motif/`.
12. `run_extract_sequences_for_phylogenetic_analyses`: The sequence for the top cooccuring motifs list are extracted for multiple sequence alignment and for phylogenetic analyses
13. `run_compress_fasta_for_phylogenetic_analyses`: Extracts phylum and kingdom information from the sequences for phylogenetic analyses and classification


