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
   
2. If not already installed, please install R and the taxize library on your local machine:
```
sudo apt-get install r-base
sudo Rscript -e 'install.packages("taxize")'
```
   
3. Install all Python dependencies (we recommend to first setup a virtual Python environment): 
```
pip3 install -r requirements.txt
```
    

# How to download the data

There are two options on how to download the data:
1. Using a data dump (includes all data from the DOOR Database and NCBI Genbank to reproduce the results from the paper)
2. Download all the data from scratch from DOOR3 and NCBI Genbank

## Download data dump
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
   
## [Alternative]: Download data from DOOR3 and NCBI Genbank

This is an alternative on how to download the data. We recommend to use the data dump for reproducibility. Downloading the data from DOOR3 and NCBI Genbank will fetch the latest data and might also include additional data that was not present during the primary analysis.

TODO PLEASE DESCRIBE



# How to run the pipeline

To run the full pipeline you just have to run the main bash script:

```
./run.sh
```

TODO: PLEASE DESCRIBE THE STEPS WHICH ARE EXECUTED BY THE BASH SCRIPT RUN AND GIVE DETAILS ABOUT THE SCRIPTS

*Individual scripts and files required to run the workflow can be found in the `bin` folder.*



