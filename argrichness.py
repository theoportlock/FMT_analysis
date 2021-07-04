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
file_card='../downstream_data/hs_10_4_igc2.CARD.tsv'
file_gene_info='../downstream_data/IGC2.1990MSPs.tsv'

# Load data
sample_type = pd.read_csv(file_sample_type, index_col=0)

gut_taxonomy = pd.read_csv(file_gut_taxonomy, index_col=0)
gut_msp_data = dd.read_csv(file_gut_msp_data)
gut_names = pd.read_table(file_gut_names, index_col=0)
card = pd.read_csv(file_card, sep='\t')
card["gene_name"] = card.ORF_ID.str.split(expand=True)[0]
card["gene_name"] = card["gene_name"].str[0:-2]
gene_info = dd.read_csv(file_gene_info, sep='\t')

sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))

gut_names["User ID"] = gut_names["User ID"].str.replace("Baseline", "Day 0")
gut_formatted_names = gut_names['User ID'].str.split(expand=True)
gut_formatted_names.loc[gut_formatted_names[0]=="Donor", 2] = gut_formatted_names.loc[gut_formatted_names[0]=="Donor", 1]
gut_formatted_names.iloc[:,2] = gut_formatted_names.iloc[:,2].str.replace("Stool","0")
gut_formatted_names = pd.merge(gut_formatted_names, sample_type, left_on=0, right_on="Patient", how='left')

replacedbygenename = gut_msp_data.merge(gene_info, how='inner', left_on='Unnamed: 0', right_on='gene_id')
replacedbycard = replacedbygenename.merge(card, how='inner', left_on='gene_name', right_on='gene_name')
del replacedbygenename

fmtcols = [col for col in replacedbycard.columns.values if 'FMT' in col]
fmtcols.append('Best_Hit_ARO')
sumofsm = replacedbycard[fmtcols].groupby('Best_Hit_ARO').sum().compute()

del replacedbycard
del gene_info
del card

gc.collect()

sumofsm.columns = pd.MultiIndex.from_frame(
        gut_formatted_names[[0,2,'Type']],
        names=['Patient','Days after treatment', 'Type'])

gut_count_df = pd.DataFrame(sumofsm.astype(bool).sum(axis=0))
gut_count_df = pd.DataFrame(sumofsm.sum(axis=0))

gut_patient_count_df = gut_count_df.drop(["Donor"],axis=0, level=0).reset_index()
gut_patient_count_df = gut_patient_count_df.rename(columns={0:"ARG richness"})
gut_patient_count_df["Days after treatment"] = gut_patient_count_df["Days after treatment"].astype(int)

gut_donor_count_df = gut_count_df.xs(["Donor"],axis=0, level=0).reset_index()
#gut_donor_count_df = pd.DataFrame(gut_msp_donor_data.astype(bool).sum(axis=0).compute())
gut_donor_count_df = gut_donor_count_df.rename(columns={0:"ARG richness"})
gut_donor_count_df.Type="Donor stool"

# Plot
x='Days after treatment'
y='ARG richness'
hue='Type'

#fig, (gut_donor_ax, gut_patient_ax, oral_patient_ax) = plt.subplots(1, 3, sharey=True, gridspec_kw={'width_ratios': [1, 7, 7]})
fig, (gut_donor_ax, gut_patient_ax) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 7]})
#fig, gut_patient_ax = plt.subplots()

sns.boxplot(ax = gut_donor_ax, data=gut_donor_count_df, x="Type", y=y, linewidth=1).set_title("donor")
sns.stripplot(ax = gut_donor_ax, data=gut_donor_count_df, x="Type", y=y, linewidth=1)
sns.boxplot(ax = gut_patient_ax, data=gut_patient_count_df, x=x, y=y, hue=hue, linewidth=1).set_title("gut")
sns.stripplot(ax = gut_patient_ax, data=gut_patient_count_df, x=x, y=y, hue=hue, linewidth=1)

# Format Plot
gut_patient_ax.legend([],[], frameon=False)
gut_donor_ax.legend([],[], frameon=False)
gut_patient_ax.set_ylabel('')
gut_donor_ax.set_xlabel('')
sns.despine(trim=True, left=True)
box_pairs=[((0,'FMT'),(7,'FMT'))]

# Stats
ax, test_results = add_stat_annotation(
    gut_patient_ax,
    data=gut_patient_count_df,
    x=x,
    y=y,
    hue=hue,
    test='Mann-Whitney',
    text_format='full',
    #loc='outside',
    box_pairs=box_pairs,
    verbose=2)

# Print
plt.tight_layout()
#plt.show()
plt.savefig("argsum.jpg")
