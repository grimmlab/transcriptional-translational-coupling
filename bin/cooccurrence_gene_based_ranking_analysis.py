# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 18:14:12 2020

@author: ge67mef
"""

import pandas as pd
import numpy as np
import os
import argparse


parser = argparse.ArgumentParser(description='co-occurrence based gene ranking')                                                              
                                                                                                                                                          
parser.add_argument('final_output_file', type=str,                                                                                                        
                    help='a final combined input file containing only both (translational and transcriptional) classified genomes')                       
                                                                                                                                                          
#parser.add_argument('path_plots' , type=str,                                                                                                              
#                    help='path for all plots generated')                                                                                                  
                                                                                                                                                          
parser.add_argument('path_outputfile', type=str,                                                                                                         
                    help='path for all output file')                                                                                                                                
args = parser.parse_args()                                                                                                                
final_output_file = args.final_output_file                                                                                                                     
path_outputfile = args.path_outputfile                           

#final_output_file = 'Final_combined_files.txt'

final_output = pd.read_csv(final_output_file, sep='\t')

def count_to_dict(it, dic):
    if it == 'na' or it == 'NA' or it == 'hypothetical protein' or it == 'conserved hypothetical protein' or it == '-':
        return
    
    if it in dic.keys():
        dic[it] += 1
    else:
        dic[it] = 1
        
    return

def count_to_dict_2(it, dic):    
    if it in dic.keys():
        dic[it] += 1
    else:
        dic[it] = 1
        
    return


def add_to_dict(it, dic, val):
    if it == 'na' or it == 'NA' or it == 'hypothetical protein' or it == 'conserved hypothetical protein' or it == '-':
        return
    
    if it in dic.keys():
        dic[it] = dic[it] + ',' + val.rstrip()
    else:
        dic[it] = val.rstrip()
    
    return


def count_val_to_dict(source_dic, target_dic, split_ch):
    
    for k in source_dic.keys():
        target_dic[k] = len(source_dic[k].split(split_ch))
        
def find_all_in_str(k,s):
    l_k = len(k.split(','))
    c = 0
    for i in k.split(','):
        if i in s:
            c = c+1
    if c == l_k:
        val = True
    else:
        val = False

    return val        

func_dict = dict()
func_dict_name = dict()
func_dict_name_cnt = dict()

gene_dict = dict()
gene_dict_name = dict()
gene_dict_name_cnt = dict()

cognum_dict = dict()
cognum_dict_name = dict()
cognum_dict_name_cnt = dict()

for i in range(0,len(final_output)):
    func_lst = final_output.iloc[i, 9].split(';')[0:-1] # function
    gene_lst = final_output.iloc[i, 7].split(',')# gene
    cognum_lst = final_output.iloc[i, 8].split(',') # cog number
    
    for itm in func_lst:
        count_to_dict(itm, func_dict)
        add_to_dict(itm, func_dict_name, final_output.iloc[i,1])
    
    for itm in gene_lst:
        count_to_dict(itm.lower(), gene_dict)
        add_to_dict(itm.lower(), gene_dict_name, final_output.iloc[i,1])    
        
    for itm in cognum_lst:
        count_to_dict(itm, cognum_dict)
        add_to_dict(itm, cognum_dict_name, final_output.iloc[i,1])    
    

count_val_to_dict(func_dict_name, func_dict_name_cnt, ',')
count_val_to_dict(gene_dict_name, gene_dict_name_cnt, ',')
count_val_to_dict(cognum_dict_name, cognum_dict_name_cnt, ',')

func_ana_output = pd.DataFrame(
    {#'start': operon_start,
      #'stop': operon_stop,
      'genome':list(func_dict_name.keys()),
      'name': list(func_dict_name.values()),
      'count': list(func_dict_name_cnt.values()),
      'genome_count': list(func_dict.values())
    })
    
#func_ana_output.to_csv('functional_occurrence_based_ranking_output.txt', sep='\t')
    
gene_ana_output = pd.DataFrame(
    {#'start': operon_start,
      #'stop': operon_stop,
      'genome':list(gene_dict_name.keys()),
      'name': list(gene_dict_name.values()),
      'count': list(gene_dict_name_cnt.values()),
      'genome_count': list(gene_dict.values())
    })    
    
#gene_ana_output.to_csv('gene_occurrence_based_ranking_output.txt', sep='\t')    
    
cognum_ana_output = pd.DataFrame(
    {#'start': operon_start,
      #'stop': operon_stop,
      'genome':list(cognum_dict_name.keys()),
      'name': list(cognum_dict_name.values()),
      'count': list(cognum_dict_name_cnt.values()),
      'genome_count': list(cognum_dict.values())
    })    
    
#cognum_ana_output.to_csv('cog_occurrence_based_ranking_output.txt', sep='\t')


# generate sorted column of gene_symbol
il=[]
for i in final_output['corrected_gene_symbol']: 
    tmp = i.rstrip().split(',')
    tmp.sort()
    il.append(','.join(tmp))
    
final_output['sorted_gene_symbol'] = il
    
final_output['sorted_gene_symbol'] = il    
gene_ana_output_sort = gene_ana_output.sort_values(ascending=False, by=['count'])
gsum = np.sum(gene_ana_output_sort['count'])

top_cutoff_num = 10;
top_gene_ana = gene_ana_output_sort[0:top_cutoff_num]
        
sorted_gsym_dict = dict()

for g in top_gene_ana['genome']:
    for i in range(0,len(final_output)):
        if g in final_output['sorted_gene_symbol'][i].lower():
            sorted_gsym_dict[i] = final_output['sorted_gene_symbol'][i].lower()


master_gsymb_list = list(sorted_gsym_dict.values())

sorted_top_gene = list(top_gene_ana['genome'])
sorted_top_gene.sort()
#sorted_top_gene = sorted_top_gene[0:4]

motif_count = dict()

add_val = dict()
motif_find = dict()
add_indices = dict()
add_s_list = dict()

sorted_gsym_l = list(final_output['sorted_gene_symbol'])
gsym_l = list(final_output['corrected_gene_symbol'])

        
        
for i in master_gsymb_list:
    il = i.split(',')
    #il=['prfA','rplQ']
    
    for s in sorted_gsym_l:
        motif = []
        for m in il:
            #t_str = ' ' + m.lower() + ' '
            #m_str = ' ' + s.lower().replace(',',' ') + ' '
            t_str = m.lower()
            m_str = s.lower().replace(',',' ')
            if t_str in m_str:
                if m not in motif:
                    motif.append(m)
                #add_val.append(sorted_gsym_l.index(s))
        if len(motif) > 1:
            motif_find[','.join(motif)] = 1
    
for k in motif_find.keys():
    for s in sorted_gsym_l:
        val = find_all_in_str(k.lower(),s.lower())
        if val:
            count_to_dict_2(k,motif_count)
            add_val[k] = sorted_gsym_l.index(s)
            add_s_list[k] = s
            
for k in add_s_list.keys():
    indices = [i for i, x in enumerate(sorted_gsym_l) if x == add_s_list[k]]
    add_indices[k] = indices           
        
#        if len(motif) > 1:
#            if not ','.join(motif) in motif_count.keys():
#                count_to_dict_2(','.join(motif),motif_count)
                

add_val_k = list(add_val.values())
coccur_stat = pd.DataFrame(
    {'motif': list(motif_count.keys()),
      'count': list(motif_count.values())  
    })    


gene_output_text =  os.path.join(args.path_outputfile, 'gene_cooccurrence_based_ranking_output.txt')
coccur_stat.to_csv(gene_output_text,index=False,sep='\t')
