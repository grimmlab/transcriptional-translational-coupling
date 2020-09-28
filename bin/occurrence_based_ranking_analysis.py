# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 18:14:12 2020

@author: ge67mef
"""

import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import os


parser = argparse.ArgumentParser(description='occurrence based ranking and visualization. ')

parser.add_argument('final_output_file', type=str,
                    help='a final combined input file containing only both (translational and transcriptional) classified genomes')

parser.add_argument('path_plots' , type=str,
                    help='path for all plots generated')

parser.add_argument('path_outputfiles', type=str,
                    help='path for all output files and plots')

args = parser.parse_args()

final_output_file = args.final_output_file
path_plots = args.path_plots
path_outputfiles = args.path_outputfiles


final_output = pd.read_csv(final_output_file, sep='\t')

def count_to_dict(it, dic):
    if it == 'na' or it == 'NA' or it == 'hypothetical protein' or it == 'conserved hypothetical protein' or it == '-':
        return
    
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
      #'name': list(func_dict_name.values()),
      'count': list(func_dict_name_cnt.values()),
      'genome_count': list(func_dict.values())
    })
   
func_output_text =  os.path.join(args.path_outputfiles, 'functional_occurrence_based_ranking_output.txt') 
func_ana_output.to_csv(func_output_text, sep='\t')
func_ana_output_sort = func_ana_output.sort_values(ascending=False, by=['count'])
gsum = np.sum(func_ana_output_sort['count'])

top_cutoff_num = 20;

top_func_ana = func_ana_output_sort[0:top_cutoff_num]
perc_func_ana = list(top_func_ana['count']/gsum * 100)
labels_fun = top_func_ana['genome']
sizes_fun = perc_func_ana


fig1, ax1 = plt.subplots()
fig1.set_size_inches(10, 10)
clrs = sns.color_palette('husl', n_colors=top_cutoff_num) 
ax1.pie(sizes_fun, labels=labels_fun, autopct='%1.1f%%', startangle=90,colors=clrs)
centre_circle = plt.Circle((0,0),0.70,fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)
ax1.tick_params(labelsize=10)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
func_output_fig1 =  os.path.join(args.path_plots, 'occurrence_functional_pie_chart.png')
fig1.savefig(func_output_fig1, dpi=300,bbox_inches='tight')
plt.show()

fig, ax = plt.subplots()
fig.set_size_inches(10, 10)
xax = range(0,len(top_func_ana))
plt.bar(xax, perc_func_ana)
xax = range(0,len(top_func_ana))
plt.xticks(xax, top_func_ana['genome'])

plt.xlabel("Funtions",fontsize=14)
plt.ylabel("Percentage occurrence",fontsize=14)
plt.tick_params(labelsize=8,rotation=90)
func_output_fig2 =  os.path.join(args.path_plots, 'occurrence_functional_bar_graph.png')
fig.savefig(func_output_fig2,dpi=300,bbox_inches='tight')
plt.show()



gene_ana_output = pd.DataFrame(
    {#'start': operon_start,
      #'stop': operon_stop,
      'genome':list(gene_dict_name.keys()),
      #'name': list(gene_dict_name.values()),
      'count': list(gene_dict_name_cnt.values()),
      'genome_count': list(gene_dict.values())
    })    
    

gene_output_text =  os.path.join(args.path_outputfiles, 'gene_occurrence_based_ranking_output.txt')
gene_ana_output.to_csv(gene_output_text, sep='\t')    
    
gene_ana_output_sort = gene_ana_output.sort_values(ascending=False, by=['count'])
gsum = np.sum(gene_ana_output_sort['count'])

top_cutoff_num = 20;

top_gene_ana = gene_ana_output_sort[0:top_cutoff_num]
#print (top_gene_ana['count'])
perc_gene_ana = list(top_gene_ana['count']/gsum * 100)
labels_gene = top_gene_ana['genome']
sizes_gene = perc_gene_ana
#print (sizes)

fig1, ax1 = plt.subplots()
fig1.set_size_inches(10, 10)
clrs = sns.color_palette('husl', n_colors=top_cutoff_num) 
ax1.pie(sizes_gene, labels=labels_gene, autopct='%1.1f%%', startangle=90,colors=clrs)
centre_circle = plt.Circle((0,0),0.70,fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)
ax1.tick_params(labelsize=10)
ax1.axis('equal')
gene_output_fig1 =  os.path.join(args.path_plots, 'occurrence_gene_pie_chart.png')
fig1.savefig(gene_output_fig1,dpi=300)
plt.show()

fig, ax = plt.subplots()
fig.set_size_inches(10, 10)
xax = range(0,len(top_gene_ana))
plt.bar(xax, perc_gene_ana)
xax = range(0,len(top_gene_ana))
plt.xticks(xax, top_gene_ana['genome'])
plt.xlabel("Genes",fontsize=14)
plt.ylabel("Percentage occurrence",fontsize=14)
plt.tick_params(labelsize=12,rotation=90)
gene_output_fig2 =  os.path.join(args.path_plots, 'occurrence_gene_bar_graph.png')
fig.savefig(gene_output_fig2,dpi=300)
plt.show()


cognum_ana_output = pd.DataFrame(
    {#'start': operon_start,
      #'stop': operon_stop,
      'genome':list(cognum_dict_name.keys()),
      #'name': list(cognum_dict_name.values()),
      'count': list(cognum_dict_name_cnt.values()),
      'genome_count': list(cognum_dict.values())
    })    
cog_output_text =  os.path.join(args.path_outputfiles, 'cog_occurrence_based_ranking_output.txt')
cognum_ana_output.to_csv(cog_output_text, sep='\t')    
cognum_ana_output_sort = cognum_ana_output.sort_values(ascending=False, by=['count'])
gsum = np.sum(cognum_ana_output_sort['count'])

top_cutoff_num = 20;

top_cog_ana = cognum_ana_output_sort[0:top_cutoff_num]
perc_cog_ana = list(top_cog_ana['count']/gsum * 100)
labels_cog = top_cog_ana['genome']
sizes_cog = perc_cog_ana


fig1, ax1 = plt.subplots()
fig1.set_size_inches(10, 10)
clrs = sns.color_palette('husl', n_colors=top_cutoff_num) 
ax1.pie(sizes_cog, labels=labels_cog, autopct='%1.1f%%', startangle=90,colors=clrs)
centre_circle = plt.Circle((0,0),0.70,fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)
ax1.tick_params(labelsize=10)
ax1.axis('equal')
cog_output_fig1 =  os.path.join(args.path_plots, 'occurrence_cog_pie_chart.png')
fig1.savefig(cog_output_fig1,dpi=300,bbox_inches='tight')
plt.show()


fig, ax = plt.subplots()
fig.set_size_inches(10, 10)
xax = range(0,len(top_cog_ana))
plt.bar(xax, perc_cog_ana)
xax = range(0,len(top_cog_ana))
plt.xticks(xax, top_cog_ana['genome'])
plt.xlabel("COG ID",fontsize=14)
plt.ylabel("Percentage occurrence",fontsize=14)
plt.tick_params(labelsize=12,rotation=90)
cog_output_fig2 =  os.path.join(args.path_plots, 'occurrence_cog_bar_graph.png')
fig.savefig(cog_output_fig2,dpi=300,bbox_inches='tight')
plt.show()

