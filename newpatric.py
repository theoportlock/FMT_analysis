#!/usr/bin/env python import gc
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import dask.dataframe as dd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

gut_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
gut_msp_data = dd.read_csv("../dowDownstream/ norm.csv")
patric = dd.read_csv("patricNorm.csv")
gene_info = dd.read_csv(file_gene_info, sep='\t')
patric_info = pd.read_csv("../downstream_data/PATRIC_sp_gene.txt", sep='\t', index_col="PATRIC ID")


replacedbygenename = gut_msp_data.merge(gene_info, how='inner', left_on='Unnamed: 0', right_on='gene_id')
replacedbypatric = replacedbygenename.merge(patric, how='inner', left_on='gene_name', right_on='igc2_id')
del replacedbygenename

fmtcols = [col for col in replacedbypatric.columns.values if 'FMT' in col]
fmtcols.append('vf_id')
sumofsm = replacedbypatric[fmtcols].groupby('vf_id').sum().compute()

del replacedbypatric, gene_info, patric

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
plt.savefig('patric.jpg')
plt.show()
