# -*- coding: utf-8 -*-
#!/usr/bin/python3

__author__ = "Richa Bharti, Dominik Grimm"
__copyright__ = "Copyright 2019-2022"
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = "Richa Bharti, Dominik Grimm"
__email__ = "richabharti74@gmail.com"
__status__ = "Dev"


import requests
import json
import pickle
import urllib
from bs4 import BeautifulSoup
from ftplib import FTP
from pathlib import Path
import gzip
import shutil
import glob
import os
import argparse
import shutil
import random



parser = argparse.ArgumentParser(description='Fetch door2 and genbank data. ')

parser.add_argument('output_path', type=str,
                    help='path where fetched data should to be stored.')

args = parser.parse_args()

output_path = args.output_path
# output_path = './'

# f = open('door2_response.pckl', 'wb')
# pickle.dump(obj, f)
# f.close()

f = open('door2_response.pckl', 'rb')
response = pickle.load(f)
f.close()

#door2_root_url = 'http://161.117.81.224/DOOR3/'
door2_ncid_url = 'http://161.117.81.224/DOOR3/organisms_ajax.php?mode=DataTable&_=1590143213326'

#response = requests.post(door2_ncid_url)

# Handling response
json_obj = json.loads(response.text)

ncids_list = json_obj['aaData']

# Collect all relevant
ncid_dict = dict()
for ncid_item in ncids_list:
    if ncid_item[0] in ncid_dict.keys():
        #print(ncid_item[0] + ' already present in the dict!')
        continue
    ncid_dict[ncid_item[0].rstrip()] = ncid_item[1]


# Filter the dict to get one (C) entry from each ncid list
ncid_filter_dict = dict()
for key in ncid_dict.keys():
    href_links = ncid_dict[key]

    # collect the links
    href_links = BeautifulSoup(href_links, "html.parser")
    href_list = []
    for link in href_links.findAll('a'):
        st = link.get('href').find('id=') + 3
        en = link.get('href').find('&')
        href_list.append(link.get('href')[st:en])

    for i in range(0,len(ncid_dict[key].split('</a>'))-1):
       t_ncid = ncid_dict[key].split('</a>')[i].split('>')[1]
       if '(C)' in t_ncid:
           ncid_filter_dict[key] =  href_list[i] + ',' + t_ncid
           break

# Save to list for future use
n_key_list = list(ncid_filter_dict.keys())
k_offset = 0
#n_key_list = n_key_list[k_offset:]

n_key_list = random.choices(n_key_list, k = 150)

# Download the operon files
for key in n_key_list:
    t_ncid = ncid_filter_dict[key].split(',')[0]
    t_nc_name = ncid_filter_dict[key].split(',')[1][0:-3]
    sp_name = key.rstrip().replace(' ','_').replace('/','').replace(':','')
    Path(output_path + '/' + sp_name).mkdir(parents=True, exist_ok=True)
    dwn_filename = sp_name + '_' + t_nc_name+'.txt'
   # dwn_url = door2_root_url+'downloadNCoperon.php?NC_id='+t_ncid
    #print(dwn_url+': '+dwn_filename)
    #print('Downloading ' + str(n_key_list.index(key)) + ': ' + sp_name + '/' + dwn_filename)
    print('Downloading ' + str(n_key_list.index(key)) + ': ' + dwn_filename)
    #urllib.request.urlretrieve(dwn_url, output_path + '/' + sp_name +'/' + dwn_filename)


# NCBI part
#assembly_summary = pd.read_csv('assembly_summary_genbank.txt', sep='\t', comment="#  ")
#assembly_list = list(assembly_summary['assembly_accession'])
#ftp_list = list(assembly_summary['ftp_path'])

f = open('assembly_accession.pckl', 'rb')
assembly_list = pickle.load(f)
f.close()

f = open('ftp_path.pckl', 'rb')
ftp_list = pickle.load(f)
f.close()

f = open('assembly_gcf_accession.pckl', 'rb')
assembly_gcf_list = pickle.load(f)
f.close()

#skip_folder = []

#skip_folder = ['Bacteriovorax marinus SJ', 'Borrelia afzelii PKo', 'Blattabacterium sp. (Blattella germanica) str. Bge', 'Aggregatibacter actinomycetemcomitans D11S-1', 'Neisseria gonorrhoeae TCDC-NG08107', 'Streptococcus salivarius 57.I', 'Borrelia turicatae 91E135', 'Chlamydophila pneumoniae AR39', 'Enterococcus hirae ATCC 9790', 'Tistrella mobilis KA081020-065', 'Shigella sonnei 53G', 'Pseudomonas fluorescens A506', 'Clostridium difficile BI1', 'Pseudomonas putida ND6', 'Escherichia coli NA114', 'Chlamydia psittaci NJ1', 'Pyrococcus abyssi GE5', 'Clostridium perfringens SM101', 'Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344', 'Comamonas testosteroni CNB-1', 'Clostridium difficile 630', 'Borrelia hermsii DAH'] 

skip_folder = ["Acinetobacter baumannii AB307-0294","Acinetobacter baumannii MDR-TJ","Aggregatibacter actinomycetemcomitans D11S-1","Arthrobacter sp. Rue61a","Bacteriovorax marinus SJ","Blattabacterium sp. (Blattella germanica) str. Bge","Borrelia afzelii PKo","Borrelia burgdorferi B31","Borrelia hermsii DAH","Borrelia turicatae 91E135","Burkholderia sp. CCGE1001","Chlamydia psittaci CP3","Chlamydia psittaci NJ1","Chlamydophila pneumoniae AR39","Clostridium difficile 630","Clostridium difficile BI1","Clostridium perfringens SM101","Comamonas testosteroni CNB-1","Enterococcus faecalis 62","Enterococcus hirae ATCC 9790","Escherichia coli NA114","Escherichia coli UMNK88","Klebsiella oxytoca E718","Lactobacillus amylovorus GRL 1112","Lactobacillus buchneri CD034","Lactobacillus plantarum subsp. plantarum ST-III","Neisseria gonorrhoeae TCDC-NG08107","Pseudomonas fluorescens A506","Pseudomonas putida ND6","Pyrococcus abyssi GE5","Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344","Shigella sonnei 53G","Streptococcus salivarius 57.I","Thermoanaerobacterium saccharolyticum JW/SL-YS485","Tistrella mobilis KA081020-065","Xanthomonas albilineans GPE PC73"]

ncbi_root_url = 'https://www.ncbi.nlm.nih.gov./assembly/?term='
ncbi_ftp_root_url = 'ftp.ncbi.nlm.nih.gov'

for key in n_key_list:
    #continue

    #print (skip_folder)
    # get the local directory name
    sp_name = key.rstrip().replace(' ','_').replace('/','').replace(':','')

    # get NC_xxx
    t_nc_name = ncid_filter_dict[key].split(',')[1][0:-3]
    #print(ncbi_root_url+t_nc_name)

    # get the GCA_xxx from the url of the response
    ncbi_response = requests.post(ncbi_root_url+t_nc_name)
    s = requests.session()
    s.keep_alive = False
    #ncbi_response.close()
    #ncbi_response.connection.close()
    #print(ncbi_response.url)

    if not 'GCF_' in ncbi_response.url:
        skip_folder.append(key)
        print('GCF not found in ncbi url, skipping: ' + key)
        continue;


    # get the GCA_xxx from the url of the response
    ftp_id = 'GCA_' + ncbi_response.url.split('/')[-2].split('_')[1]
    # ftp_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/'+ ftp_id[0:3] + '/' + ftp_id[3:6] + '/' + ftp_id[6:9]

    # get the ftp directory path from the ftp_path column list
    if ftp_id in assembly_list:
        ftp_ind = assembly_list.index(ftp_id)
    elif 'GCF' + ftp_id[3:] in assembly_gcf_list:
          print('Replacing GCA key with GCF key: ' + ftp_id)
          ftp_ind = assembly_gcf_list.index('GCF' + ftp_id[3:])
    else:
          skip_folder.append(key)
          print('ftp id :' + ftp_id + ' not found, skipping: ' + key)
          #ftp_ls.close()
          continue;

    #print(ftp_id + ' ' + str(n_key_list.index(key)+k_offset) + ': ' + ftp_list[ftp_ind])

    # establish ftp connection to the ncbi ftp server
    ftp_ls = FTP(ncbi_ftp_root_url)

    # perform an anonymous login the ncbi ftp server
    ftp_ls.login()

    #print (ftp_id, str(ftp_ind), ftp_list[ftp_ind])

    st = ftp_list[ftp_ind].find('gov/') + 4
    t_dir = ftp_list[ftp_ind][st:]
    #print(t_dir)

    # navigate the corresponding GCA genomes directory
    ftp_ls.cwd(t_dir)

    # fetch the exact file name which is to be downloaded
    ftp_file_list = ftp_ls.nlst()
    ftp_targ_file = []
    for s in ftp_file_list:
        if 'feature_table' in s or 'genomic.gff.gz' in s or 'genomic.gtf.gz' in s or 'protein.faa.gz' in s or  'genomic.fna.gz' in s:
            ftp_targ_file.append(s)
    if len(ftp_targ_file) == 0:
          skip_folder.append(key)
          print('None of the required files found in the ftp server, skipping: ' + key)
          ftp_ls.close()
          continue;

    #ftp_cnt_file = ftp_ls.nlst('*feature_table*')
    ftp_cnt_file = ftp_targ_file

    # open a local file with the same file name
    #print (ftp_cnt_file[0])
    
    for dl_file in ftp_cnt_file:
    
        loc_file = open(output_path + '/' + sp_name + '/' + dl_file, 'wb')
        #loc_file = open('tmp.txt.gz', 'wb')

        # perform a ftp get on the file which is to be downloaded in chunks of 1024 bytes
        print("Downloading " + str(n_key_list.index(key)+k_offset) + ' ' + dl_file + ' from the ncbi ftp server...')
        ftp_ls.retrbinary('RETR ' + dl_file, loc_file.write, 1024)
        
        # close the downloaded file to save it
        loc_file.close()

    # close the ftp connection (otherwise it becomes unstable and timesout)
    ftp_ls.quit()
    ftp_ls.close()

    # del ftp_ls
    
    


# unzip and delete zip file
for key in n_key_list:

    # get the local directory name
    sp_name = key.rstrip().replace(' ','_').replace('/','').replace(':','')

    if key in skip_folder:
        print ('Deleting... ' + sp_name)
        if os.path.isdir(output_path + '/' + sp_name):
            shutil.rmtree(output_path + '/' + sp_name)        
        continue

    # get the gz file name
    gz_file = glob.glob( output_path + '/' + sp_name+"/*.gz")
        
    for gz_f in gz_file:

        # get the txt file name
        txt_file = gz_f[0:-3]   
        print ('Gunzip ' + str(n_key_list.index(key)) + ':' + gz_f)

        # unzipping operation (gunzip)
        with gzip.open(gz_f, 'rb') as f_in:
            with open(txt_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

            f_out.close()

        # remove the gz file
        print('Removing ' + gz_f)
        os.remove(gz_f)
        
        
# Run prokka 
for key in n_key_list:

    # get the local directory name
    sp_name = key.rstrip().replace(' ','_').replace('/','').replace(':','')
    
    if key in skip_folder:
        continue    

    # get the gz file name
    prokka_input_protien = glob.glob( output_path + '/' + sp_name+"/*_protein.faa")
    prokka_input_genomic_ls = glob.glob( output_path + '/' + sp_name+"/*_genomic.fna")
    
    for p_in in prokka_input_genomic_ls:
        if 'rna_from' in p_in or 'cds_from' in p_in:
            continue
        else:
            prokka_input_genomic = p_in
            
    prokka_input_protien = prokka_input_protien[0]
    prokka_output_folder = output_path + '/' + sp_name + '/' + 'prokka_output'
    cmd = 'prokka --kingdom Bacteria --proteins ' + prokka_input_protien + ' ' + prokka_input_genomic + ' --outdir ' + prokka_output_folder + ' --usegenus --compliant --cpus 64 --rfam'
    print(cmd)
    os.system(cmd)
