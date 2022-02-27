# -*- coding: utf-8 -*-
#!/usr/bin/python3

__author__ = "Richa Bharti"
__copyright__ = "Copyright 2019"
__credits__ = ["Richa Bharti"]
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = "Richa Bharti, Dominik Grimm"
__email__ = "richabharti74@gmail.com"
__status__ = "Dev"

import pandas as pd
import numpy as np
import argparse
import os


parser = argparse.ArgumentParser(description='co-occurrence based functional ranking ')                                                      
                                                                                                                                                  
parser.add_argument('final_output_file', type=str,                                                                                                
                    help='a final combined input file containing only both (translational and transcriptional) classified genomes')               
                                                                                                                                                  
#parser.add_argument('path_plots' , type=str,                                                                                                     
#                    help='path for all plots generated')                                                                                         
                                                                                                                                                  
parser.add_argument('path_outputfile', type=str,                                                                                                 
                    help='path for all output file')                                                                                   
                                                                                                                                                  
args = parser.parse_args()                                                                                                                        
                                                                                                                                                  
final_output_file = args.final_output_file                                                                                                        
#path_plots = args.path_plots                                                                                                                     
path_outputfile = args.path_outputfile                                                                                                          
                                                                                                                                                  
#final_output_file = 'Final_combined_files.txt'                                                                                                   
                                                                                                                                                  

final_output = pd.read_csv(final_output_file, sep='\t')
final_output['corrected_function'] = final_output['corrected_function'].str.lower()

def count_to_dict(it, dic):
    if it =='na' or it == 'NA' or it == 'hypothetical protein' or it == 'conserved hypothetical protein' or it == '-':
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
    if it =='na' or it == 'NA' or it == 'hypothetical protein' or it == 'conserved hypothetical protein' or it == '-':
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
    l_k = len(k.split(';'))
    c = 0
    for i in k.split(';'):
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
  
il=[]
for i in final_output['corrected_function']: 
    tmp = i.rstrip().split(';')[1:-1]
    tmp.sort()
    il.append(';'.join(tmp))
  
final_output['sorted_func_symbol'] = il    


gene_ana_output_sort = gene_ana_output.sort_values(ascending=False, by=['count'])
gsum = np.sum(gene_ana_output_sort['count'])

func_ana_output_sort = func_ana_output.sort_values(ascending=False, by=['count'])
fsum = np.sum(func_ana_output_sort['count'])



top_cutoff_num = 20;
top_func_ana = func_ana_output_sort[0:top_cutoff_num]
        
sorted_func_dict = dict()

for g in top_func_ana['genome']:
    for i in range(0,len(final_output)):
        if g in final_output['sorted_func_symbol'][i]:
            sorted_func_dict[i] = final_output['sorted_func_symbol'][i]


master_func_list = list(sorted_func_dict.values())

sorted_top_func = list(top_func_ana['genome'])
sorted_top_func.sort()
#sorted_top_gene = sorted_top_gene[0:4]

motif_count = dict()
motif_genename_ind = dict()

add_val = dict()
motif_find = dict()
add_indices = dict()
add_s_list = dict()

#final_output['sorted_func_symbol'] = final_output['sorted_func_symbol'].str.lower()


sorted_func_l = list(final_output['sorted_func_symbol'])
func_l = list(final_output['corrected_function'])

        
        
for i in master_func_list:
    il = i.split(';')
    #il=['prfA','rplQ']
    
    for s in sorted_func_l:
        motif = []
        for m in il:
            #t_str = ' ' + m.lower() + ' '
            #m_str = ' ' + s.lower().replace(';',' ') + ' '
            t_str = m.lower()
            m_str = s.lower().replace(';',' ')
            if t_str in m_str:
                if m not in motif:
                    motif.append(m)
                #add_val.append(sorted_gsym_l.index(s))
        if len(motif) > 1:
            motif_find[';'.join(motif)] = 1
    


for k in motif_find.keys():
    ind=0
    for s in sorted_func_l:
        val = find_all_in_str(k.lower(),s.lower())
        #ind = sorted_func_l.index(s)
        if val:
            count_to_dict_2(k,motif_count)
            add_val[k] = sorted_func_l.index(s)
            add_s_list[k] = s
            add_to_dict(k,motif_genename_ind,str(ind))
        ind = ind + 1
            
            
for k in add_s_list.keys():
    indices = [i for i, x in enumerate(sorted_func_l) if x == add_s_list[k]]
    add_indices[k] = indices           
        
#        if len(motif) > 1:
#            if not ','.join(motif) in motif_count.keys():
#                count_to_dict_2(','.join(motif),motif_count)
                

add_val_k = list(add_val.values())


pre_coccur_stat = pd.DataFrame(
    {'motif': list(motif_count.keys()),
      'count': list(motif_count.values()),
      'indices':list(motif_genename_ind.values())
    })
    
pre_coccur_stat.sort_values(ascending=False, inplace = True, by=['count'])


top_counts = 20
motif_genename = []
cnt=1
for indl in pre_coccur_stat['indices'][0:top_counts]:
    t_str = []
    ncid_str = []
    for ind in indl.split(','):
        gstr = final_output['Filename '][int(ind)].rstrip()
        t_str.append(gstr.split('NC')[0].replace('_',' ').rstrip())
        ncid_str.append('NC'+gstr.split('NC')[1])
    motif_genename.append(','.join(t_str))
    genome_list = None
    genome_list = pd.DataFrame(
            {'genome':t_str,
             'ncid':ncid_str
            })
    top_genomes_file = os.path.join(args.path_outputfile+'/top_genome_list/'+'Rank_' + str(cnt) + str('_genomelist.txt'))
    print(top_genomes_file)
    genome_list.to_csv(top_genomes_file,sep='\t',index=False)
    cnt = cnt + 1


#coccur_genename_stat = pd.DataFrame()
#coccur_genename_stat['count'] = pre_coccur_stat['count'][0:top_counts-1]
#coccur_genename_stat['motif'] = pre_coccur_stat['motif'][0:top_counts-1]
#coccur_genename_stat['genome'] = list(motif_genename)
#coccur_genename_stat.to_csv('cooccurrence_func_statistics_motif_with_genomename.txt', sep='\t')

coccur_stat = pd.DataFrame(
    {'motif': list(motif_count.keys()),
      'count': list(motif_count.values())  
    })    

                                                                                                                
func_output_text =  os.path.join(args.path_outputfile, 'functional_cooccurrence_based_ranking_output.txt')     
coccur_stat.to_csv(func_output_text,index=False,sep='\t')                                                       
