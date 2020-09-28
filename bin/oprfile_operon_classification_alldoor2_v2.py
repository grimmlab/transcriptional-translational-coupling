import textwrap
import pandas as pd
import numpy as np
import json
import argparse



parser = argparse.ArgumentParser(description='Operon classification. ')

parser.add_argument('final_tab_file', type=str,
                    help='a feature table file from genbank')

parser.add_argument('opr_file', type=str,
                    help='a operon file from door2 database')

args = parser.parse_args()

fin_tab_file = args.final_tab_file
opr_file = args.opr_file

#fin_tab_file = 'Acetobacter_pasteurianus_IFO_3283-01-42C\GCA_000010945.1_ASM1094v1_feature_table.txt'
#opr_file = 'Acetobacter_pasteurianus_IFO_3283-01-42C\Acetobacter_pasteurianus_IFO_3283-01-42C_NC_017150.txt'

#fin_tab_file = 'test/GCA_000005845.2_ASM584v2_feature_table.txt'
#opr_file = 'test/Escherichia_coli strK-12_substr_MG1655_NC_000913(C).opr'


print('Feature file: ' + fin_tab_file.split('/')[-1])
print('Operon file: ' + opr_file.split('/')[-1])

def op_classify_gene(gene_name):
    if (gene_name in ref_rib_genes):
        gene_type = 'translational'
    elif (gene_name in ref_rnapolym_genes):
        gene_type = 'transcriptional'
    else: 
        gene_type = 'none'

    return gene_type

def op_classify_name(name):
    if any([s for s in ref_rib_func if s.lower() in name]) or (any([s for s in name if s in ref_rib_genes])):
        gene_type = 'translational'
    elif any([s for s in ref_rnapolym_func if s.lower() in name]) or (any([s for s in name if s in ref_rnapolym_genes])):
        gene_type = 'transcriptional'
    else:
        gene_type = 'none'

    return gene_type

nc_file = pd.read_csv(opr_file, sep='\t', comment='<')
 
nc_dict = dict()
 
for i in range(0,len(nc_file)):
    nc_dict[nc_file['Synonym'][i]] = nc_file['Product'][i]


locus_gene_dict = dict()
locus_name_dict = dict()

final_tab = pd.read_csv(fin_tab_file, sep='\t')

if final_tab['locus_tag'].isnull().values.all():
#    #print (opr_file.split('/')[-1])
    print ('\t Cannot process! Missing locus_tag column!' )
    # print (opr_file.split('/')[-2].replace('_',' ')  + '\t' + opr_file.split('/')[-1])
    exit()

final_tab['symbol'] = final_tab['symbol'].replace(np.nan, 'no_symbol', regex=True)
final_tab['name'] = final_tab['name'].replace(np.nan,'hypothetical protein', regex=True)

for i in range(0,len(final_tab)):
   
   if 'gene' == final_tab['# feature'][i]:
       continue
   locus_gene_dict[final_tab['locus_tag'][i]] = final_tab['symbol'][i]
   if final_tab['locus_tag'][i] in nc_dict.keys():
       locus_name_dict[final_tab['locus_tag'][i]] = nc_dict[final_tab['locus_tag'][i]]
   else:
       locus_name_dict[final_tab['locus_tag'][i]] = 'hypothetical protein'

opr_input = pd.read_csv(opr_file,sep='\t', comment='<')
opr_input['Synonym'] = opr_input[['OperonID','Synonym']].groupby(['OperonID'])['Synonym'].transform(lambda x: ','.join(x))
opr_input['COG_number'] = opr_input[['OperonID','COG_number']].groupby(['OperonID'])['COG_number'].transform(lambda x: ','.join(x))
opr_input = opr_input.drop_duplicates(subset='OperonID', keep='first') # Remove duplicate rows

operon_start = []
operon_stop = []
operon_list = []
operon_cognum = []
for l in opr_input.index:
   operon_list.append(opr_input['Synonym'][l])
   operon_start.append(opr_input['Start'][l])
   operon_stop.append(opr_input['End'][l])
   operon_cognum.append(opr_input['COG_number'][l])

rnapolym_input = pd.read_csv('RNApolymerase_F.txt', sep='\t')
ref_rnapolym_genes = [x.lower() for x in rnapolym_input["Gene"]]

ribosomes_input = pd.read_csv('Ribosomes_newlist_F.txt', sep='\t')
ref_rib_genes = [x.lower() for x in ribosomes_input["Gene"]]

ref_rib_func = ['Ribosomal RNA small subunit','Ribosomal RNA large subunit',
                  'Translation initiation factor','Ribosome-binding factor',
                  '30S', '50S','23S','5S','Prolyl-tRNA','Cysteinyl-tRNA','Leucyl-tRNA',
                  'Glutaminyl-tRNA','Seryl-tRNA','Asparaginyl-tRNA',
                  'Phenylalanyl-tRNA','16S','Valyl-tRNA', 'Arginyl-tRNA',
                  'Aspartyl-tRNA','Peptidyl-tRNA','Alanyl-tRNA',
                  'Glutamyl-tRNA','Histidyl-tRNA','Anticodon:']


ref_rnapolym_func = ['RNA polymerase', 'transcriptional regulator',
                'Transcription', 'RNAP subunit beta',
                'RNA polymerase sigma factor', 'DNA-directed RNA polymerase subunit',
                'sigma54 factor', 'Sigma-70',' Sigma-38','sigma']

operon_class = []
prot_class = []
gene_name_list = []
for i in range(0,len(operon_list)):
#for i in range(52,53):
    o_list = operon_list[i].lstrip()
    o_list = o_list.split(',')   
    o_len = len(o_list)

    op_type_list = []
    op_name_list = []
    g_name_l = ''
    for j in range(0,o_len):
        if o_list[j] in locus_gene_dict.keys():
            if 'no_symbol' == locus_gene_dict[o_list[j]]:
                #print ('no gene name found')
                if o_list[j] in locus_name_dict.keys():
                    op_type = op_classify_name(locus_name_dict[o_list[j]].lower())
                else:
                    op_type = 'none'  
                op_gene = 'NA'
            else:
                #print (locus_gene_dict[o_list[j]] + ' found')
                op_type = op_classify_gene(locus_gene_dict[o_list[j]].lower())
                if 'none' == op_type:
                    if o_list[j] in locus_name_dict.keys():
                        op_type = op_classify_name(locus_name_dict[o_list[j]].lower())
                    else:
                        op_type = 'none'              
                op_gene = locus_gene_dict[o_list[j]]
        else:
            #print (o_list[j] + ' key not found')
            op_type = 'none'
            op_gene = 'NA'

        op_type_list.append(op_type)
        
        if o_list[j] in locus_name_dict.keys():
            op_name_list.append(locus_name_dict[o_list[j]] + ';')
        else:
            op_name_list.append('none')

        g_name_l = g_name_l  + op_gene + ','

    n_translational = op_type_list.count('translational') 
    n_transcriptional = op_type_list.count('transcriptional')

    if (n_translational > 0) and (n_transcriptional > 0):
        operon_class.append('both')
    elif (n_translational > 0):
        operon_class.append('translation')
    elif (n_transcriptional > 0):
        operon_class.append('transcription')
    else:
        operon_class.append('none')

    prot_class.append(''.join(op_name_list))

    gene_name_list.append(g_name_l[0:-1])

final_operon_list = pd.DataFrame(
    {#'start': operon_start,
      #'stop': operon_stop,
      'function':prot_class,
      'operon': operon_list,
      'classification': operon_class,
      'gene_symbol': gene_name_list,
      'COG_number': operon_cognum
    })

output_file = opr_file[0:-4] + '_operon_classification_output.txt'
final_operon_list.to_csv(output_file, sep='\t')

