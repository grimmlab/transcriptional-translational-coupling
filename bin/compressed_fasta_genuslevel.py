# -*- coding: utf-8 -*-
#!/usr/bin/python3

__author__ = “Richa Bharti”
__copyright__ = “Copyright 2019”
__credits__ = [“Richa Bharti”]
__license__ = “MIT”
__version__ = “0.1.0”
__maintainer__ = “Richa Bharti, Dominik Grimm”
__email__ = “richabharti74@gmail.com”
__status__ = “Dev”

from Bio import SeqIO
import argparse
import re
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
txz = importr('taxize')


parser = argparse.ArgumentParser(description='extracted uncompressed 16S rRNA fasta file')                                                  
                                                                                                                                         
parser.add_argument('fasta_file', type=str,                                                                                       
                    help='extracted ncbi 16S rRNA fasta file for a motif')      
                                                                                                                                         
parser.add_argument('path_outputfile', type=str,
                    help='path for all output file')

args = parser.parse_args()                                                                                                               
                                                                                                                                         
input_file = args.fasta_file                                                                                               
path_outputfile = args.path_outputfile
#input_file = 'nusA_infB_genome_sequences.fasta'
output_file = args.path_outputfile+'/'+input_file.split('/')[-1].split('.')[0] + '_reduced.fasta'

w_handle = open(output_file,'w')


records = list(SeqIO.parse(input_file, "fasta"))

id_dict = dict()
for i in range(0,len(records)):
    wrd = records[i].description.split(' ')[1]
    
    if wrd in id_dict.keys():
        continue
    else:
        id_dict[wrd] = 1
        nid = records[i].name 
        #id_line = '>' + ' '.join(records[i].description.split(' ')[0:2]) + '\n'
        genus_id = ' '.join(records[i].description.split(' ')[1:2])
        # Call R function
        tr = txz.tax_name(sci=genus_id, get = ["phylum", "kingdom"], db = "itis", messages = False)
        id_line = '>' + tr[2][0] + '|' + tr[3][0] + '|' + tr[1][0]  + '\n'
        seq = records[i]._seq._data + '\n'
        w_handle.write(id_line)
        seq_nl =  re.sub("(.{80})", "\\1\n", seq, 0, re.DOTALL)
        w_handle.write(seq_nl)

w_handle.close()
        
        
