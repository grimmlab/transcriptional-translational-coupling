
# Reannotation check pipeline

This commands will fetch 150 genomes for reannotation with prokka. We will then compare these reannotations with those from GenBank.

1. Install Prokka
```
conda install -c conda-forge -c bioconda -c defaults prokka
```

3. Uncompress the ftp_path to ftp_path.pckl file before running the analyses
4. Run the analyses with the command below 

```
unzip ftp_path.zip
mkdir results
python3 fetch_150_random_genbank_ids_restapi_run_prokka.py results
```
