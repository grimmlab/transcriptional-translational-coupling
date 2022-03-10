# -*- coding: utf-8 -*-
#!/usr/bin/python3.6

__author__ = "Richa Bharti"
__copyright__ = "Copyright 2019"
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = "Richa Bharti, Dominik Grimm"
__email__ = "richabharti74@gmail.com"
__status__ = "Dev"

import textwrap
import pandas as pd
import numpy as np
import json
import argparse
import matplotlib.pyplot as plt
import itertools 
import collections
import os
#nusA, infB
#rpoA, rplQ, rpsD, rpsM, rpsK
#rplA, rplK
parser = argparse.ArgumentParser(description='gene motif search')                                                      
                                                                                                                                                  
parser.add_argument('final_input_file', type=str,                                                                                                
                    help='a final combined input file containing only both (translational and transcriptional) classified genomes')               
                                                                                                                                                  
#parser.add_argument('path_plots' , type=str,                                                                                                     
#                    help='path for all plots generated')                                                                                         
                                                                                                                                                  
parser.add_argument('path_outputfile', type=str,                                                                                                 
                    help='path for all output file')                                                                                   
                                                                                                                                                  
parser.add_argument('input_gene_motif', type=str,                                                                                                 
                    help='input gene combination eg. nusG, nusA, rpoB, secY')                                                                                   

args = parser.parse_args()                                                                                                                        


input_gene = args.input_gene_motif #'nusG, nusA, rpoB, secY'
final_output_file = args.final_input_file #'Final_combined_files_corrected.txt'

final_output = pd.read_csv(final_output_file, sep='\t')

#final_output[' function '] = final_output[' function '].str.lower()


def count_to_dict(it, dic):    
    if it in dic.keys():
        dic[it] += 1
    else:
        dic[it] = 1
        
    return


def add_to_dict(key, dic, val):
    
    if key in dic.keys():
        dic[key] = dic[key] + ',' + val.rstrip()
    else:
        dic[key] = val.rstrip()
    
    return


gene_symbol_sorted = []

for i in final_output['corrected_gene_symbol']: 
    tmp = i.rstrip().split(',')
    tmp.sort()
    gene_symbol_sorted.append(','.join(tmp).lower())
    
gene_dict = dict()
genome_dict = dict()
ind = 0    
for i in gene_symbol_sorted:
    count_to_dict(i, gene_dict)
    add_to_dict(i, genome_dict, final_output['Filename '][ind])
    ind = ind + 1
    
    
input_gene_sorted = input_gene.replace(' ', '').rstrip().lower().split(',')
input_gene_key = ','.join(input_gene_sorted).lower()
gl = len(input_gene_sorted)


gsum = 0
gfound = []
genome_found = []
gfound_count = []
for k in gene_dict.keys():    
    gc = 0
    for g in input_gene_sorted:
        if g in k:
            gc = gc + 1
            continue
        else:
            break
        
    if gc == gl:
        #print(k, gene_dict[k])
        gsum = gsum + gene_dict[k]
        gfound.append(k)
        gfound_count.append(gene_dict[k])
        genome_found.append(genome_dict[k])
        

final_genome_ls = []
final_ncid_ls = []
final_gfound = []
final_gfound_count = []
c=0
for i in genome_found:
    t_str = []
    ncid_str = []
    gf_tmp = []
    gf_c_tmp = []
    gf_tmp.append(gfound[c])
    gf_c_tmp.append(gfound_count[c])            
    for g in i.split(','):
        t_str.append(g.split('NC')[0].replace('_',' ').rstrip())
        ncid_str.append('NC'+g.split('NC')[1])
        gf_tmp.append('')
        gf_c_tmp.append('')        
    final_genome_ls.extend(t_str)
    final_ncid_ls.extend(ncid_str)
    final_gfound.extend(gf_tmp[0:-1])
    final_gfound_count.extend(gf_c_tmp[0:-1])
    c = c + 1
        
print('Input gene: ' + input_gene + ' is found ' + str(gsum) + ' times!')

genome_csv = pd.DataFrame(
    {'genome': final_genome_ls,
     'ncid': final_ncid_ls,
     'motif': final_gfound,
     'count':final_gfound_count
    })    
filename =  os.path.join(args.path_outputfile, input_gene.replace(' ','').replace(',','_') +'_genome_list.txt')     
#filename = input_gene.replace(' ','').replace(',','_') + '_genome_list.txt'
genome_csv.to_csv(filename, sep='\t', index=False)

genome_found_csv = pd.DataFrame(
    {'motif': gfound,
     'count': gfound_count
    })  
filename_count =  os.path.join(args.path_outputfile, input_gene.replace(' ','').replace(',','_') +'_motif_list.txt')     
#filename_count = input_gene.replace(' ','').replace(',','_') + '_motif_list.txt'
genome_found_csv.to_csv(filename_count, sep='\t', index=False)
# if input_gene_key in gene_dict.keys():
#     print(gene_dict[input_gene_key])
# else:
#     print('Gene sybmol:' + input_gene + ' cannot be found!')
