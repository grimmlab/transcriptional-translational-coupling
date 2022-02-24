![Python 3.7](https://img.shields.io/badge/Python-3.7-blue.svg) ![Ubuntu 18.04.2](https://img.shields.io/badge/Ubuntu-18.04.2-green.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Systematic analysis of the underlying genomic architecture for transcriptional-translational coupling in prokaryotes


This is the associated GitHub page to reproduce the results from the following paper (under preparation):


**Systematic analysis of the underlying genomic architecture for transcriptional-translational coupling in prokaryotes**  
Richa Bharti, Daniel Siebert, Bastian Blombach, Dominik G. Grimm
   

## Pipeline Summary
We created a comparative genomics pipeline for screening genomic distributions of probable conserved operonic motifs in bacteria (Figure 1). More details can be found in the accompanying manuscript.

<p align="center">
  <img src="https://github.com/grimmlab/transcriptional-translational-coupling/blob/master/Figure%201.png">
</p>


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


