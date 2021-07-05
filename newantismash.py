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
file_antismash='../downstream_data/igc2.antismash.simple.txt'
file_gene_info='../downstream_data/IGC2.1990MSPs.tsv'

# Load data
sample_type = pd.read_csv(file_sample_type, index_col=0)

gut_taxonomy = pd.read_csv(file_gut_taxonomy, index_col=0)
gut_msp_data = dd.read_csv(file_gut_msp_data)
gut_names = pd.read_table(file_gut_names, index_col=0)
antismash = pd.read_table(file_antismash)
antismash.columns = ['gene_name','sm','id']
gene_info = dd.read_csv(file_gene_info, sep='\t')

sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))

gut_names["User ID"] = gut_names["User ID"].str.replace("Baseline", "Day 0")
gut_formatted_names = gut_names['User ID'].str.split(expand=True)
gut_formatted_names.loc[gut_formatted_names[0]=="Donor", 2] = gut_formatted_names.loc[gut_formatted_names[0]=="Donor", 1]
gut_formatted_names.iloc[:,2] = gut_formatted_names.iloc[:,2].str.replace("Stool","0")
gut_formatted_names = pd.merge(gut_formatted_names, sample_type, left_on=0, right_on="Patient", how='left')

replacedbygenename = gut_msp_data.merge(gene_info, how='inner', left_on='Unnamed: 0', right_on='gene_id')
replacedbyantismash = replacedbygenename.merge(antismash, how='inner', left_on='gene_name', right_on='gene_name')
del replacedbygenename

fmtcols = [col for col in replacedbyantismash.columns.values if 'FMT' in col]
fmtcols.append('sm')
sumofsm = replacedbyantismash[fmtcols].groupby('sm').sum().compute()

del replacedbyantismash, gene_info, antismash

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
plt.legend(title='sm', bbox_to_anchor=(1.001, 1), loc='upper left', fontsize='small')

# gut fmt
att1 = sumofsm.drop('Donor',level=0,axis='columns')
att3 = att1.xs('FMT',level=2,axis='columns').groupby(level=1,axis = 1).mean()
att3['0']=att3['0'].div(att3['0'].sum(axis=0), axis=0)
att3['7']=att3['7'].div(att3['7'].sum(axis=0), axis=0)
att3['30']=att3['30'].div(att3['30'].sum(axis=0), axis=0)
att3['90']=att3['90'].div(att3['90'].sum(axis=0), axis=0)
att3.columns = att3.columns.astype(int)
att3 = att3.reindex(columns=att3.columns.sort_values())
sns.heatmap(att3,cmap='binary',ax=gut_patient_fmt_ax,cbar=False)
#att3.T.plot(kind='bar',stacked=True, ax=gut_patient_fmt_ax, legend=False)
gut_patient_fmt_ax.title.set_text('FMT - Stool')

# gut placebo
att1 = sumofsm.drop('Donor',level=0,axis='columns')
att3 = att1.xs('PLACEBO',level=2,axis='columns').groupby(level=1,axis = 1).mean()
att3['0']=att3['0'].div(att3['0'].sum(axis=0), axis=0)
att3['7']=att3['7'].div(att3['7'].sum(axis=0), axis=0)
att3['30']=att3['30'].div(att3['30'].sum(axis=0), axis=0)
att3['90']=att3['90'].div(att3['90'].sum(axis=0), axis=0)
att3.columns = att3.columns.astype(int)
att3 = att3.reindex(columns=att3.columns.sort_values())
sns.heatmap(att3,cmap='binary',ax=gut_patient_placebo_ax)
#att3.T.plot(kind='bar',stacked=True, ax=gut_patient_placebo_ax)
gut_patient_placebo_ax.title.set_text('Placebo - Stool')

plt.setp(gut_patient_placebo_ax.get_legend().get_texts(), fontsize = 8)
plt.setp(gut_patient_placebo_ax.get_legend().get_title(), fontsize = 8)
plt.tight_layout()
plt.savefig('secmetab.jpg')
#plt.show()
