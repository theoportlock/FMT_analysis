#!/usr/bin/env python
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

# Data location
file_gut_taxonomy="../downstream_data/taxo.csv"
file_oral_taxonomy="../oral_downstream_data/oraltaxo.csv"
file_gut_msp_data="../downstream_data/merged.final.mgs.med.vec.10M.csv"
#file_oral_msp_data="../oral_downstream_data/merged.final.mgs.med.vec.10M.csv"
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

# Replace column headers with multiindex
sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))

gut_names["User ID"] = gut_names["User ID"].str.replace("Baseline", "Day 0")
gut_formatted_names = gut_names['User ID'].str.split(expand=True)
gut_formatted_names.loc[gut_formatted_names[0]=="Donor", 2] = gut_formatted_names.loc[gut_formatted_names[0]=="Donor", 1]
gut_formatted_names.iloc[:,2] = gut_formatted_names.iloc[:,2].str.replace("Stool","0")
gut_formatted_names = pd.merge(gut_formatted_names, sample_type, left_on=0, right_on="Patient", how='left')
gut_formatted_names[2] = gut_formatted_names[2].astype(int)
gut_msp_data.columns = pd.MultiIndex.from_frame(
        gut_formatted_names[[0,2,'Type']],
        names=['Patient','Days after treatment', 'Type'])

oral_names["User ID"] = oral_names["User ID"].str.replace("Baseline", "Day 0")
oral_formatted_names = oral_names['User ID'].str.split(expand=True)
oral_formatted_names.loc[oral_formatted_names[0]=="Donor", 2] = oral_formatted_names.loc[oral_formatted_names[0]=="Donor", 1]
oral_formatted_names.iloc[:,2] = oral_formatted_names.iloc[:,2].str.replace("Stool","0")
oral_formatted_names = pd.merge(oral_formatted_names, sample_type, left_on=0, right_on="Patient", how='left')
oral_formatted_names[2] = oral_formatted_names[2].astype(int)
oral_msp_data.columns = pd.MultiIndex.from_frame(
        oral_formatted_names[[0,2,'Type']],
        names=['Patient','Days after treatment', 'Type'])

# Count number of MSPs
gut_msp_patient_data = gut_msp_data.drop("Donor", axis=1, level=0)
gut_msp_patient_data = gut_msp_patient_data.drop("PLACEBO", axis=1, level=2)
gut_msp_patient_data.sort_index(axis=1,level=1,inplace=True)
gut_msp_patient_data.columns = gut_msp_patient_data.columns.to_flat_index()

gut_msp_donor_data = gut_msp_data.xs("Donor", axis=1, level=0)
gut_msp_donor_data.columns = gut_msp_donor_data.columns.droplevel(1)

oral_msp_patient_data = oral_msp_data.drop("Donor", axis=1, level=0)
oral_msp_patient_data = oral_msp_patient_data.drop("PLACEBO", axis=1, level=2)
oral_msp_patient_data.sort_index(axis=1,level=1,inplace=True)
oral_msp_patient_data.columns = oral_msp_patient_data.columns.to_flat_index()

# Plot
fig, (gut_donor_ax, gut_patient_ax, oral_patient_ax) = plt.subplots(1, 3, sharey=True, gridspec_kw={'width_ratios': [1, 7, 7]})

# Format Plot
gut_patient_ax.legend([],[], frameon=False)
gut_donor_ax.legend([],[], frameon=False)
gut_patient_ax.set_ylabel('')
oral_patient_ax.set_ylabel('')
gut_donor_ax.set_xlabel('')
sns.despine(trim=True, left=True)
box_pairs=[((0,'FMT'),(7,'FMT'))]

# Print
plt.tight_layout()
sns.heatmap(gut_msp_patient_data.T,cmap='viridis',vmax=0.00005,yticklabels=True,xticklabels=False)
plt.show()
#plt.savefig("results/species_richness.pdf")
