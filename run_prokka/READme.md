
# Reannotation check pipeline

This commands will fetch 150 genomes for reannotation with prokka. We will then compare these reannotations with those from GenBank.

1. Install Prokka
```
conda install -c conda-forge -c bioconda -c defaults prokka
```

Alternatively you can install Prokka as follows:  
```
sudo apt-get install libdatetime-perl libxml-simple-perl libdigest-perl-md5-perl git default-jre bioperl
sudo cpan Bio::Perl
git clone https://github.com/tseemann/prokka.git 
prokka/bin/prokka --setupdb
```

More detailed information can be found here: https://github.com/tseemann/prokka



2. Uncompress the ftp_path to ftp_path.pckl file before running the analyses  

```
unzip ftp_path.zip
```

4. Run the analyses with the command below 

```
mkdir results
python3 fetch_150_random_genbank_ids_restapi_run_prokka.py results
```
