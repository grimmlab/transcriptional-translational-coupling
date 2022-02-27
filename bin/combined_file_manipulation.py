__author__ = “Richa Bharti”
__copyright__ = “Copyright 2019”
__credits__ = [“Richa Bharti”]
__license__ = “MIT”
__version__ = “0.1.0”
__maintainer__ = “Richa Bharti, Dominik Grimm”
__email__ = “richabharti74@gmail.com”
__status__ = “Dev”


import textwrap
import pandas as pd
import numpy as np
import json
import argparse


parser = argparse.ArgumentParser(description='Modifying classification file. ')                                                    
                                                                                                                                                
parser.add_argument('final_output_file', type=str,                                                                                              
                    help='a final combined input file containing only both (translational and transcriptional) classified genomes')             
                                                                                                                                                
                                                                                                                                                
args = parser.parse_args()                                                                                                                      
                                                                                                                                                
fin_tab_file = args.final_output_file                                                                                                      

final_tab = pd.read_csv(fin_tab_file, sep='\t')

gene_cog_dict = dict()
cog_gene_dict = dict()
for i in range(0,len(final_tab)):
    s = final_tab[' gene_symbol '][i]
    c = final_tab[' COG_number'][i]
    
    gs = s.strip().split(',')
    cn = c.strip().split(',')
   
    if len(gs) == len(cn):
        for ig in range(0,len(gs)):
            if gs[ig].strip() == 'NA' or cn[ig].strip() == '-':
                continue
            else:
                gene_cog_dict[gs[ig].strip()] = cn[ig].strip()
                cog_gene_dict[cn[ig].strip()] = gs[ig].strip()
    else:
        print (str(i))


g_final = []
c_final = []
for i in range(0,len(final_tab)):        
    s = final_tab[' gene_symbol '][i]
    c = final_tab[' COG_number'][i]
    
    gs = s.strip().split(',')
    cn = c.strip().split(',')
   
    g_ls = []
    c_ls = []
    if len(gs) == len(cn):
        for ig in range(0,len(gs)):
            if gs[ig].strip() == 'NA':
                if cn[ig] in cog_gene_dict.keys():
                    g_ls.append(cog_gene_dict[cn[ig].strip()])
                    c_ls.append(cn[ig].strip())
                else:
                    g_ls.append(gs[ig].strip())
                    c_ls.append(cn[ig].strip())
                continue
            
            if cn[ig].strip() == '-':
                if gs[ig] in gene_cog_dict.keys():
                    c_ls.append(gene_cog_dict[gs[ig].strip()])
                    g_ls.append(gs[ig].strip())
                else:
                    c_ls.append(cn[ig].strip())
                    g_ls.append(gs[ig].strip())
                continue
            g_ls.append(gs[ig].strip())
            c_ls.append(cn[ig].strip())
    else:
        print (str(i))
        for g in gs:
            g_ls.append(g.strip())
            
        for cgn in cn:
            c_ls.append(cgn.strip())
            
    g_final.append(','.join(g_ls))
    c_final.append(','.join(c_ls))



gene_fn_dict = dict()    
for i in range(0,len(final_tab)):
    s = g_final[i]
    f = final_tab[' function '][i]
    
    gs = s.strip().split(',')
    fn = f.strip().split(';')[:-1]
    
    if len(gs) == len(fn):
        for ig in range(0,len(gs)):
            if gs[ig].strip() != 'NA':
                gene_fn_dict[gs[ig].strip()] = fn[ig].strip()
    else:
        print (str(i))
                
f_final = []        
for i in range(0,len(final_tab)):
    s = g_final[i]
    f = final_tab[' function '][i]
    
    gs = s.strip().split(',')
    fn = f.strip().split(';')[:-1]
    
    f_ls = []    
    if len(gs) == len(fn):
        for ig in range(0,len(gs)):
            if gs[ig].strip() != 'NA':
                f_ls.append(gene_fn_dict[gs[ig].strip()])
            else:
                f_ls.append(fn[ig].strip())
                
                
    else:
        for fnc in fn:
            f_ls.append(fnc.strip())        
            
    f_final.append(';'.join(f_ls) + ';')        
            
    
final_tab['corrected_gene_symbol'] = g_final
final_tab['corrected_cog_nummber'] = c_final    
final_tab['corrected_function'] = f_final


final_tab.to_csv(fin_tab_file[:-4]+'_corrected.txt', sep='\t')
