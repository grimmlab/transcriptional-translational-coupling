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
    

# How to download the data

We created a data dump, including all necessary data from the DOOR3 database and NCBI Genbank to reproduce the results from the paper. Alternatively all data can be also fetched from the DOOR3 database and from NCBI Genbank (time consuming).

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

In the following we will give some detailes about the individual steps of the pipeline within the `run.sh` bash script:

1. The environment and all PATH variables are setup by the script.
2. `run_uncompress_door2_and_ncbi_data`: The data dumps from the GitHub repo are merged and unzipped
3. `run_classification_code`: This function is running the 


