![Python 3.7](https://img.shields.io/badge/Python-3.7-blue.svg) ![Ubuntu 18.04.2](https://img.shields.io/badge/Ubuntu-18.04.2-green.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Systematic analysis of the underlying genomic architecture for transcriptional-translational coupling in prokaryotes


**Systematic analysis of the underlying genomic architecture for transcriptional-translational coupling in prokaryotes**  
Richa Bharti, Daniel Siebert, Bastian Blombach, Dominik G. Grimm

 <p style='text-align: justify;'> The aim of this study is to provide a comprehensive workflow to systematically investigate bacterial genomes for the abundance of transcriptional and translational associated genes clustered in distinct operons.</p>

 ## Pipeline Summary
We have created a comparative genomics pipeline for screening genomic distributions of probable conserved operonic motifs in bacteria (Figure 1). More details can be found in the accompanying manuscript. The workflow provided here contain series of steps written as invidual python scripts which are stiched together and can be controlled by bash scripts `run.sh`.

This is the associated GitHub page can be used to reproduce the results in tablular and some figures from the following paper (under preparation):


<p align="center">
  <img src="https://github.com/grimmlab/transcriptional-translational-coupling/blob/master/Figure%201.png">
</p>

## Prerequisites
Basic prerequisites that need to be satisfied:

### OS
Any Linux based distro should work. We tested the scripts using:

Distributor ID: Ubuntu <br/>
Description:    Ubuntu 18.04.2 LTS <br/>
Release:        18.04 <br/>
Codename:       bionic <br/>

`lsb_release -a` on a Ubuntu based system.

###  Dependencies and packages
<p style='text-align: justify;'> In order to reproduce the results please make sure that you have all the libaries on your local machine.
Without these the workflow will fail to run. </p>

- R has to be installed before you install all the Python dependencies `sudo apt-get install r-base`   

- Install all Python dependencies (we recommend to first setup a virtual Python environment): `pip3 install -r requirements.txt`   

### How to download the data
Data will be automatically fetched by the accompanion bash script when running the pipeline. 
Downliading the data from the DOOR2 database and GenBank can take several hours.

TODO: PLEASE DESCRIBE THIS IN MORE DETAIL

### How to run the pipeline

To run the full pipeline you just have to run the main bash script:

```
./run.sh
```
 
TODO: PLEASE DESCRIBE THE STEPS WHICH ARE EXECUTED BY THE BASH SCRIPT RUN AND GIVE DETAILS ABOUT THE SCRIPTS

*Individual scripts and files required to run the workflow can be found in the `bin` folder.*



