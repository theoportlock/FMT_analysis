#!/usr/bin/env python
import gc
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import dask.dataframe as dd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

# Data location
file_gut_taxonomy="../downstream_data/taxo.csv"
file_gut_msp_data="../project/Downstream/norm.csv"
file_gut_names="../downstream_data/P15952_20_03_sample_info_stool.txt"
file_sample_type="../downstream_data/PROFITplaceboTab.csv"
file_kegg='../downstream_data/hs_10_4_igc2_annot_20180614_vs_kegg82.best.csv'
file_gene_info='../downstream_data/IGC2.1990MSPs.tsv'

# Load data
sample_type = pd.read_csv(file_sample_type, index_col=0)

gut_taxonomy = pd.read_csv(file_gut_taxonomy, index_col=0)
gut_msp_data = dd.read_csv(file_gut_msp_data)
gut_names = pd.read_table(file_gut_names, index_col=0)
kegg = pd.read_csv(file_kegg)
kegg.columns = ['gene_id','gene_name','ko', 'bitscore']
gene_info = dd.read_csv(file_gene_info, sep='\t')

sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))

gut_names["User ID"] = gut_names["User ID"].str.replace("Baseline", "Day 0")
gut_formatted_names = gut_names['User ID'].str.split(expand=True)
gut_formatted_names.loc[gut_formatted_names[0]=="Donor", 2] = gut_formatted_names.loc[gut_formatted_names[0]=="Donor", 1]
gut_formatted_names.iloc[:,2] = gut_formatted_names.iloc[:,2].str.replace("Stool","0")
gut_formatted_names = pd.merge(gut_formatted_names, sample_type, left_on=0, right_on="Patient", how='left')

replacedbykegg = gut_msp_data.merge(kegg, how='inner', left_on='Unnamed: 0', right_on='gene_id')
#replacedbykegg = replacedbygenename.merge(kegg, how='inner', left_on='gene_name', right_on='gene_name')

fmtcols = [col for col in replacedbykegg.columns.values if 'FMT' in col]
fmtcols.append('ko')
sumofsm = replacedbykegg[fmtcols].groupby('ko').sum().compute()

del replacedbykegg
del gene_info
del kegg

gc.collect()

sumofsm.columns = pd.MultiIndex.from_frame(
        gut_formatted_names[[0,2,'Type']],
        names=['Patient','Days after treatment', 'Type'])

# plot
fig, (gut_donor_ax, gut_patient_fmt_ax, gut_patient_placebo_ax) = plt.subplots(1, 3, sharey=True, gridspec_kw={'width_ratios': [1, 4, 4]})

# donor
datt1 = sumofsm.xs('Donor', level=0, axis='columns').copy()
datt2 = datt1.groupby(datt1.index).sum().mean(axis=1).copy()
datt3 = pd.DataFrame(datt2)
datt3[0]=datt3[0].div(datt3[0].sum(axis=0), axis=0)
#datt3.T.plot.bar(stacked=True, ax=gut_donor_ax, legend=False)
sns.heatmap(datt3,cmap='binary',ax=gut_donor_ax, cbar=False)
gut_donor_ax.title.set_text('Donor')
gut_donor_ax.set_ylabel("Relative Abundance")
#plt.legend(title='sm', bbox_to_anchor=(1.001, 1), loc='upper left', fontsize='small')

# gut fmt
fatt1 = sumofsm.drop('Donor',level=0,axis='columns').copy()
fatt3 = fatt1.xs('FMT',level=2,axis='columns').groupby(level=1,axis = 1).mean()
fatt3['0']=fatt3['0'].div(fatt3['0'].sum(axis=0), axis=0)
fatt3['7']=fatt3['7'].div(fatt3['7'].sum(axis=0), axis=0)
fatt3['30']=fatt3['30'].div(fatt3['30'].sum(axis=0), axis=0)
fatt3['90']=fatt3['90'].div(fatt3['90'].sum(axis=0), axis=0)
fatt3.columns = fatt3.columns.astype(int)
fatt3 = fatt3.reindex(columns=fatt3.columns.sort_values())
sns.heatmap(fatt3,cmap='binary',ax=gut_patient_fmt_ax,cbar=False)
#att3.T.plot(kind='bar',stacked=True, ax=gut_patient_fmt_ax, legend=False)
gut_patient_fmt_ax.title.set_text('FMT - Stool')

# gut placebo
patt1 = sumofsm.drop('Donor',level=0,axis='columns').copy()
patt3 = patt1.xs('PLACEBO',level=2,axis='columns').groupby(level=1,axis = 1).mean()
patt3['0']=patt3['0'].div(patt3['0'].sum(axis=0), axis=0)
patt3['7']=patt3['7'].div(patt3['7'].sum(axis=0), axis=0)
patt3['30']=patt3['30'].div(patt3['30'].sum(axis=0), axis=0)
patt3['90']=patt3['90'].div(patt3['90'].sum(axis=0), axis=0)
patt3.columns = patt3.columns.astype(int)
patt3 = patt3.reindex(columns=patt3.columns.sort_values())
sns.heatmap(patt3,cmap='binary',ax=gut_patient_placebo_ax)
#patt3.T.plot(kind='bar',stacked=True, ax=gut_patient_placebo_ax)
gut_patient_placebo_ax.title.set_text('Placebo - Stool')

#plt.setp(gut_patient_placebo_ax.get_legend().get_texts(), fontsize = 8)
#plt.setp(gut_patient_placebo_ax.get_legend().get_title(), fontsize = 8)
#plt.tight_layout()
#plt.savefig('secmetab.pdf')
plt.show()
