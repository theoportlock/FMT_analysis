#!/usr/bin/env python
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

import skbio
from skbio.stats.ordination import pcoa
from scipy.spatial import distance

# Data location
file_gut_taxonomy="../downstream_data/taxo.csv"
file_oral_taxonomy="../oral_downstream_data/oraltaxo.csv"
file_gut_msp_data="../downstream_data/merged.final.mgs.med.vec.10M.csv"
file_oral_msp_data="../oral_downstream_data/merged.final.mgs.med.vec.10M.csv"
file_gut_msp_data="../downstream_data/gctDown10m.csv"
file_gut_names="../downstream_data/P15952_20_03_sample_info_stool.txt"
file_oral_names="../oral_downstream_data/A.Mardinoglu_20_05_sample_info.txt"
file_sample_type="../downstream_data/PROFITplaceboTab.csv"
file_antibio="../../rscripts/anti1.csv"

# Load data
sample_type = pd.read_csv(file_sample_type, index_col=0)
gut_taxonomy = pd.read_csv(file_gut_taxonomy, index_col=0)
gut_msp_data = pd.read_csv(file_gut_msp_data, index_col=0)
gut_names = pd.read_table(file_gut_names, index_col=0)
oral_taxonomy = pd.read_csv(file_oral_taxonomy, index_col=0)
oral_msp_data = pd.read_csv(file_oral_msp_data, index_col=0)
oral_names = pd.read_table(file_oral_names, index_col=0)
gut_gene_data = dd.read_csv(file_gut_msp_data)
antibio = pd.read_csv(file_antibio, index_col=0)

gut_gene_data = gut_gene_data.set_index("Unnamed: 0")

# Replace column headers with multiindex
sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))
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
oral_msp_data = oral_msp_data.drop('Donor',level=0,axis='columns')

# Compute bray distances
fig, (gut_ax, oral_ax) = plt.subplots(1, 2, sharey=True)
gut_ax.set_title("Gut")
oral_ax.set_title("Oral")
gut_ax.legend([],[], frameon=False)
oral_ax.legend([],[], frameon=False)

df = gut_msp_data.copy()
df.columns = df.columns.to_flat_index()
df = df.T
df = df[~df.index.duplicated(keep='last')]
Ar_dist = distance.squareform(distance.pdist(df, metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
PCoA = skbio.stats.ordination.pcoa(DM_dist)
tmp = df.index.to_numpy()
patientnums = [i[0] for i in tmp]
pn = pd.DataFrame(patientnums)
ptype = pn[0].str.replace("P.*","Patient").tolist()
day = [i[1] for i in tmp]
pday = pd.DataFrame(day)
pday = pday[0].astype(int)
pday[pday > 100] = 100
pday=pday.to_list()
addedsamples = PCoA.samples.copy()
addedsamples.set_index(pd.MultiIndex.from_tuples(df.index, names=['Patient','time','type']),inplace=True)
addedsamples.drop('Donor',level=0, axis='index').reset_index().groupby("Patient").plot(x="PC1",y="PC2",ax = gut_ax,legend=False, color='gray', linewidth=0.5, label='_nolegend_')
sns.scatterplot(data=PCoA.samples, x='PC1', y='PC2', hue=pday, style=ptype, ax=gut_ax, palette='coolwarm')

df = oral_msp_data.copy()
df.columns = df.columns.to_flat_index()
df = df.T
df = df[~df.index.duplicated(keep='last')]
Ar_dist = distance.squareform(distance.pdist(df, metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
PCoA = skbio.stats.ordination.pcoa(DM_dist)
tmp = df.index.to_numpy()
day = [i[1] for i in tmp]
pday = pd.DataFrame(day)
pday = pday[0].astype(int)
pday[pday > 100] = 100
pday=pday.to_list()
addedsamples = PCoA.samples.copy()
addedsamples.set_index(pd.MultiIndex.from_tuples(df.index, names=['Patient','time','type']),inplace=True)
addedsamples.reset_index().groupby("Patient").plot(x="PC1",y="PC2",ax = oral_ax,legend=False, color='gray', linewidth=0.5, label='_nolegend_')
sns.scatterplot(data=PCoA.samples, x='PC1',y='PC2', hue=pday, ax=oral_ax, palette='coolwarm')

#plt.tight_layout()
#plt.savefig("results/PCoA.pdf")
plt.show()
