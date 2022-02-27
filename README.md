![Python 3.7](https://img.shields.io/badge/Python-3.7-blue.svg) ![Ubuntu 18.04.2](https://img.shields.io/badge/Ubuntu-18.04.2-green.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Systematic analysis of the underlying genomic architecture for transcriptional-translational coupling in prokaryotes
 The aim of this study () is to provide a comprehensive workflow to systematically investigated bacterial genomes for the abundance of transcriptional and translational associated genes clustered in distinct operons. 

**Systematic analysis of the underlying genomic architecture for transcriptional-translational coupling in prokaryotes**  
Richa Bharti, Daniel Siebert, Bastian Blombach, Dominik G. Grimm
   

## Pipeline Summary
We have created a comparative genomics pipeline for screening genomic distributions of probable conserved operonic motifs in bacteria (Figure 1). More details can be found in the accompanying manuscript. The workflow provided here contain series of steps written as invidual python scripts which are stiched together and can be controlled by bash scripts: `run.sh`. 

This is the associated GitHub page to reproduce the results from the following paper (under preparation):


<p align="center">
  <img src="https://github.com/grimmlab/transcriptional-translational-coupling/blob/master/Figure%201.png">
</p>


## How to run the pipeline

To run the full pipeline you just have to run the main bash script:

```
./run.sh
```
  
*Individual scripts and files required to run the workflow can be found in the `bin` folder.

## Operating System

Any Linux based distro should work. We tested the scripts using:  
  
Distributor ID: Ubuntu  
Description: Ubuntu 18.04.2 LTS  
Release: 18.04  
Codename: bionic  
  
'lsb_release -a' on a Ubuntu based system.  


