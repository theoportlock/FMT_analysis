#!/usr/bin/env python
import numpy as np
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
taxa_type='species'

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
#fig, (gut_donor_ax, gut_patient_fmt_ax, gut_patient_placebo_ax, oral_patient_fmt_ax, oral_patient_placebo_ax) = plt.subplots(1, 5, sharey=True, gridspec_kw={'width_ratios': [1, 4, 4, 4, 4]})

# gut fmt
taxa_not_sum = gut_msp_data.join(gut_taxonomy[taxa_type]).set_index(taxa_type)
taxa_not_sum.columns=pd.MultiIndex.from_tuples(taxa_not_sum.columns,names=['Patient','Day','Type'])
gutfmtatt1 = taxa_not_sum.drop('Donor',level=0,axis='columns')
gutfmtatt1.drop(["P0004","P0001","P0016","P0013","P0022","P0024","P0017","P0023","P0029"],axis=1, inplace=True)

# gut placebo
taxa_not_sum = gut_msp_data.join(gut_taxonomy[taxa_type]).set_index(taxa_type)
taxa_not_sum.columns=pd.MultiIndex.from_tuples(taxa_not_sum.columns,names=['Patient','Day','Type'])
gutplacatt1 = taxa_not_sum.drop('Donor',level=0,axis='columns').copy()
gutplacatt1.drop(["P0004","P0001","P0016","P0013","P0022","P0024","P0017","P0023","P0029"],axis=1, inplace=True)

# oral fmt
taxa_not_sum = oral_msp_data.join(oral_taxonomy[taxa_type]).set_index(taxa_type)
taxa_not_sum.columns=pd.MultiIndex.from_tuples(taxa_not_sum.columns,names=['Patient','Day','Type'])
oralfmtatt1 = taxa_not_sum.drop('Donor',level=0,axis='columns').iloc[:,4:].copy()
oralfmtatt1.drop(["P0004","P0001","P0016","P0013","P0022","P0024","P0017","P0023","P0029"],axis=1, inplace=True)

# oral placebo
taxa_not_sum = oral_msp_data.join(oral_taxonomy[taxa_type]).set_index(taxa_type)
taxa_not_sum.columns=pd.MultiIndex.from_tuples(taxa_not_sum.columns,names=['Patient','Day','Type'])
oralplacatt1 = taxa_not_sum.drop('Donor',level=0,axis='columns').iloc[:,4:].copy()
oralplacatt1.drop(["P0004","P0001","P0016","P0013","P0022","P0024","P0017","P0023","P0029"],axis=1, inplace=True)

shared = pd.DataFrame()
shared['total'] = gutfmtatt1.astype(bool).sum()
numshared = []
gutfmtatt1 = gutfmtatt1.T
gutplacatt1 = gutplacatt1.T
oralfmtatt1 = oralfmtatt1.T
oralplacatt1 = oralplacatt1.T

for i in oralfmtatt1.index:
    numshared.append(oralfmtatt1.loc[i][oralfmtatt1.loc[i].astype(bool)].index.isin(gutfmtatt1.loc[i][gutfmtatt1.loc[i].astype(bool)].index).sum())

shared['numshared'] = numshared
shared['PSS'] = shared['numshared'].div(shared['total'])
shared.reset_index(inplace=True)
sns.boxplot(data=shared, x='Day', y='PSS', hue='Type')

#plt.show()
#plt.tight_layout()
plt.savefig('results/shared_species.pdf')
