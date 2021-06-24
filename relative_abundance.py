#!/usr/bin/env python import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

# Data location
file_gut_taxonomy="../downstream_data/taxo.csv"
file_oral_taxonomy="../oral_downstream_data/oraltaxo.csv"
file_gut_msp_data="../downstream_data/merged.final.mgs.med.vec.10M.csv"
file_oral_msp_data="../oral_merged_downstream_data/merged.final.mgs.med.vec.10M.csv"
file_gut_names="../downstream_data/P15952_20_03_sample_info_stool.txt"
file_oral_names="../oral_downstream_data/A.Mardinoglu_20_05_sample_info.txt"
file_sample_type="../downstream_data/PROFITplaceboTab.csv"

# Load data
sample_type = pd.read_csv(file_sample_type, index_col=0)
gut_taxonomy = pd.read_csv(file_gut_taxonomy, index_col=0)
gut_msp_data = pd.read_csv(file_gut_msp_data, index_col=0)
gut_names = pd.read_table(file_gut_names, index_col=0)
oral_taxonomy = pd.read_csv(file_oral_taxonomy, index_col=0)
oral_msp_data = pd.read_csv(file_oral_msp_data, index_col=0)
oral_names = pd.read_table(file_oral_names, index_col=0)
sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))

# Format data
taxa_type='phylum'

gut_names["User ID"] = gut_names["User ID"].str.replace("Baseline", "Day 0")
gut_formatted_names = gut_names['User ID'].str.split(expand=True)
gut_formatted_names.loc[gut_formatted_names[0]=="Donor", 2] = gut_formatted_names.loc[gut_formatted_names[0]=="Donor", 1]
gut_formatted_names.iloc[:,2] = gut_formatted_names.iloc[:,2].str.replace("Stool","0")
gut_formatted_names = pd.merge(gut_formatted_names, sample_type, left_on=0, right_on="Patient", how='left')
gut_msp_data.columns = pd.MultiIndex.from_frame(
        gut_formatted_names[[0,2,'Type']],
        names=['Patient','Days after treatment', 'Type'])

oral_names["User ID"] = oral_names["User ID"].str.replace("Baseline", "Day 0")
oral_formatted_names = oral_names['User ID'].str.split(expand=True)
oral_formatted_names.loc[oral_formatted_names[0]=="Donor", 2] = oral_formatted_names.loc[oral_formatted_names[0]=="Donor", 1]
oral_formatted_names.iloc[:,2] = oral_formatted_names.iloc[:,2].str.replace("Stool","0")
oral_formatted_names = pd.merge(oral_formatted_names, sample_type, left_on=0, right_on="Patient", how='left')
oral_msp_data.columns = pd.MultiIndex.from_frame(
        oral_formatted_names[[0,2,'Type']],
        names=['Patient','Days after treatment', 'Type'])

# plot
fig, (gut_donor_ax, gut_patient_fmt_ax, gut_patient_placebo_ax, oral_patient_fmt_ax, oral_patient_placebo_ax) = plt.subplots(1, 5, sharey=True, gridspec_kw={'width_ratios': [1, 4, 4, 4, 4]})

# donor
dtaxa_not_sum = gut_msp_data.join(gut_taxonomy[taxa_type]).set_index(taxa_type)
dtaxa_not_sum.columns=pd.MultiIndex.from_tuples(dtaxa_not_sum.columns,names=['Patient','Day','Type'])
datt1 = dtaxa_not_sum.xs('Donor', level=0, axis='columns').copy()
datt2 = datt1.groupby(datt1.index).sum().mean(axis=1).copy()
datt3 = pd.DataFrame(datt2)
datt3[0]=datt3[0].div(datt3[0].sum(axis=0), axis=0)
datt3.T.plot.bar(stacked=True, ax=gut_donor_ax, legend=False)
gut_donor_ax.title.set_text('Donor')
gut_donor_ax.set_ylabel("Relative Abundance")
plt.legend(title=taxa_type, bbox_to_anchor=(1.001, 1), loc='upper left', fontsize='small')

# gut fmt
taxa_not_sum = gut_msp_data.join(gut_taxonomy[taxa_type]).set_index(taxa_type)
taxa_not_sum.columns=pd.MultiIndex.from_tuples(taxa_not_sum.columns,names=['Patient','Day','Type'])
att1 = taxa_not_sum.drop('Donor',level=0,axis='columns')
att2 = att1.xs('FMT',level=2,axis='columns').groupby(level=1,axis = 1).mean()
att3 = att2.reset_index().groupby(taxa_type).sum().copy()
att3['0']=att3['0'].div(att3['0'].sum(axis=0), axis=0)
att3['7']=att3['7'].div(att3['7'].sum(axis=0), axis=0)
att3['30']=att3['30'].div(att3['30'].sum(axis=0), axis=0)
att3['90']=att3['90'].div(att3['90'].sum(axis=0), axis=0)
att3.columns = att3.columns.astype(int)
att3 = att3.reindex(columns=att3.columns.sort_values())
att3.T.plot(kind='bar',stacked=True, ax=gut_patient_fmt_ax, legend=False)
gut_patient_fmt_ax.title.set_text('FMT - Stool')

# gut placebo
taxa_not_sum = gut_msp_data.join(gut_taxonomy[taxa_type]).set_index(taxa_type)
taxa_not_sum.columns=pd.MultiIndex.from_tuples(taxa_not_sum.columns,names=['Patient','Day','Type'])
att1 = taxa_not_sum.drop('Donor',level=0,axis='columns')
att2 = att1.xs('PLACEBO',level=2,axis='columns').groupby(level=1,axis = 1).mean()
att3 = att2.reset_index().groupby(taxa_type).sum().copy()
att3['0']=att3['0'].div(att3['0'].sum(axis=0), axis=0)
att3['7']=att3['7'].div(att3['7'].sum(axis=0), axis=0)
att3['30']=att3['30'].div(att3['30'].sum(axis=0), axis=0)
att3['90']=att3['90'].div(att3['90'].sum(axis=0), axis=0)
att3.columns = att3.columns.astype(int)
att3 = att3.reindex(columns=att3.columns.sort_values())
att3.T.plot(kind='bar',stacked=True, ax=gut_patient_placebo_ax)
gut_patient_placebo_ax.title.set_text('Placebo - Stool')

# oral fmt
taxa_not_sum = oral_msp_data.join(oral_taxonomy[taxa_type]).set_index(taxa_type)
taxa_not_sum.columns=pd.MultiIndex.from_tuples(taxa_not_sum.columns,names=['Patient','Day','Type'])
att1 = taxa_not_sum.drop('Donor',level=0,axis='columns')
att2 = att1.xs('FMT',level=2,axis='columns').groupby(level=1,axis = 1).mean()
att3 = att2.reset_index().groupby(taxa_type).sum().copy()
att3['0']=att3['0'].div(att3['0'].sum(axis=0), axis=0)
att3['7']=att3['7'].div(att3['7'].sum(axis=0), axis=0)
att3['30']=att3['30'].div(att3['30'].sum(axis=0), axis=0)
att3['90']=att3['90'].div(att3['90'].sum(axis=0), axis=0)
att3.columns = att3.columns.astype(int)
att3 = att3.reindex(columns=att3.columns.sort_values())
att3.T.plot(kind='bar',stacked=True, ax=oral_patient_fmt_ax, legend=False)
oral_patient_fmt_ax.title.set_text('FMT - Saliva')

# oral placebo
taxa_not_sum = oral_msp_data.join(oral_taxonomy[taxa_type]).set_index(taxa_type)
taxa_not_sum.columns=pd.MultiIndex.from_tuples(taxa_not_sum.columns,names=['Patient','Day','Type'])
att1 = taxa_not_sum.drop('Donor',level=0,axis='columns')
att2 = att1.xs('PLACEBO',level=2,axis='columns').groupby(level=1,axis = 1).mean()
att3 = att2.reset_index().groupby(taxa_type).sum().copy()
att3['0']=att3['0'].div(att3['0'].sum(axis=0), axis=0)
att3['7']=att3['7'].div(att3['7'].sum(axis=0), axis=0)
att3['30']=att3['30'].div(att3['30'].sum(axis=0), axis=0)
att3['90']=att3['90'].div(att3['90'].sum(axis=0), axis=0)
att3.columns = att3.columns.astype(int)
att3 = att3.reindex(columns=att3.columns.sort_values())
att3.T.plot(kind='bar',stacked=True, ax=oral_patient_placebo_ax)
oral_patient_placebo_ax.title.set_text('Placebo - Saliva')
 
plt.setp(oral_patient_placebo_ax.get_legend().get_texts(), fontsize = 8)
plt.setp(oral_patient_placebo_ax.get_legend().get_title(), fontsize = 8)
plt.setp(gut_patient_placebo_ax.get_legend().get_texts(), fontsize = 8)
plt.setp(gut_patient_placebo_ax.get_legend().get_title(), fontsize = 8)
plt.show()
#plt.tight_layout()
#plt.savefig('results/relative_abundance.pdf')
