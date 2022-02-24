![Python 3.7](https://img.shields.io/badge/Python-3.7-blue.svg) ![Ubuntu 18.04.2](https://img.shields.io/badge/Ubuntu-18.04.2-green.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Systematic analysis of the underlying genomic architecture for transcriptional-translational coupling in prokaryotes


This is the associated GitHub page to reproduce the results from the following paper (under preparation):


**Systematic analysis of the underlying genomic architecture for transcriptional-translational coupling in prokaryotes**  
Richa Bharti, Daniel Siebert, Bastian Blombach, Dominik G. Grimm

## Abstract
Transcriptional-translational coupling is accepted to be a fundamental mechanism of gene expression in prokaryotes and therefore has been analyzed in detail. However, the underlying genomic architecture of the expression machinery has not been well investigated so far. In this study, we established a bioinformatics pipeline to systematically investigated more than 1800 bacterial genomes for the abundance of transcriptional and translational associated genes clustered in distinct operons. We identified three highly frequent operonic motifs containing transcriptional and translational genes i.e., *rplk-nusG* (motif 1; in 553 genomes), *rpoA-rplQ-rpsD-rpsK-rpsM* (motif 2; in 656 genomes) and *nusA-infB* (motif 3; in 877 genomes). Interestingly, each of the three motifs harbors a gene (*nusG*, *rpsD* and *nusA*) encoding a protein which links transcription and translation in bacteria. Phylogenetic analyses suggest an enrichment of these motifs in pathogenic bacterial phyla with more than 70% for motif 3 (i.e. *Neisseria*, *Salmonella*, and *Escherichia*) and more than 50% for motif 1 (i.e. *Treponema*, *Prevotella*, *Leptospira* and *Fusobacterium*) and motif 2 (i.e. *Helicobacter*, *Campylobacter*, *Treponema* and *Prevotella*). These insights form the basis to analyze the transcriptional regulatory mechanisms orchestrating transcriptional-translational coupling and might open novel avenues for future biotechnological approaches. 

## How to run the pipeline

To run the full pipeline you just have to run the main bash script:

```
./run.sh
```
  
Individual scripts can be found in the `bin` folder.

## Operating System

Any Linux based distro should work. We tested the scripts using:  
  
Distributor ID: Ubuntu  
Description: Ubuntu 18.04.2 LTS  
Release: 18.04  
Codename: bionic  
  
'lsb_release -a' on a Ubuntu based system.  


