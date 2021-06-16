#!/usr/bin/env python
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import skbio
from statannot import add_stat_annotation
from skbio.stats.ordination import pcoa
from scipy.spatial import distance

# Data location
file_gut_taxonomy="../downstream_data/taxo.csv"
file_oral_taxonomy="../oral_downstream_data/oraltaxo.csv"
file_gut_msp_data="../downstream_data/merged.final.mgs.med.vec.10M.csv"
file_oral_msp_data="../oral_downstream_data/merged.final.mgs.med.vec.10M.csv"
file_gut_names="../downstream_data/P15952_20_03_sample_info_stool.txt"
file_oral_names="../oral_downstream_data/A.Mardinoglu_20_05_sample_info.txt"
file_sample_type="../downstream_data/PROFITplaceboTab.csv"
file_day0="../downstream_data/PROFITbaseline.csv"
file_day7="../downstream_data/PROFITday7.csv"
file_day30="../downstream_data/PROFITday30.csv"
file_day90="../downstream_data/PROFITday90.csv"
file_siafile="../downstream_data/normalised.faecal.cytokines.markers.xlsx"

# Load data
sample_type = pd.read_csv(file_sample_type, index_col=0)
gut_taxonomy = pd.read_csv(file_gut_taxonomy, index_col=0)
gut_msp_data = pd.read_csv(file_gut_msp_data, index_col=0)
gut_names = pd.read_table(file_gut_names, index_col=0)
oral_taxonomy = pd.read_csv(file_oral_taxonomy, index_col=0)
oral_msp_data = pd.read_csv(file_oral_msp_data, index_col=0)
oral_names = pd.read_table(file_oral_names, index_col=0)
day0=pd.read_csv(file_day0)
day7=pd.read_csv(file_day7)
day30=pd.read_csv(file_day30)
day90=pd.read_csv(file_day90)
siafile = pd.read_excel(file_siafile)

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

# Format metadata
day0["Days after treatment"]=0
day7["Days after treatment"]=7
day30["Days after treatment"]=30
day90["Days after treatment"]=90
vals = pd.concat([day0,day7,day30,day90])
vals["Patient_Id"] = vals.Patient.str.replace("P01","P").str.replace("\ .*", "", regex=True)
siafile["id"] = siafile.id.str.replace("P01","P").str.replace("\ .*", "", regex=True)
siafile["time_point"] = siafile.time_point.str.replace("D","")
siafile['time_point'] = siafile.time_point.astype('int64')
merged_metadata = pd.merge(siafile, vals,  how='inner', left_on=['id','time_point'], right_on = ['Patient_Id','Days after treatment'])

# Compute bray distances and plot
fig, (gut_ax, oral_ax) = plt.subplots(1, 2, sharey=True)
gut_ax.set_title("Gut")
oral_ax.set_title("Oral")

df = gut_msp_data.copy()
df.columns = df.columns.to_flat_index()
df = df.T
df = df[~df.index.duplicated(keep='last')]
Ar_dist = distance.squareform(distance.pdist(df, metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
PCoA = skbio.stats.ordination.pcoa(DM_dist)
addedsamples = PCoA.samples.copy()
addedsamples.set_index(pd.MultiIndex.from_tuples(df.index, names=['Patient','time','type']),inplace=True)
addedsamples.reset_index(inplace=True)
addedsamples.time = addedsamples.time.astype('int')
gut_merged_metadata = pd.merge(merged_metadata, addedsamples, how='inner', left_on=['id','time_point'], right_on=['Patient', 'time'])
sns.scatterplot(data=gut_merged_metadata, x='PC1', y='PC2', hue='PPI', ax=gut_ax, palette='colorblind')

df = oral_msp_data.copy()
df.columns = df.columns.to_flat_index()
df = df.T
df = df[~df.index.duplicated(keep='last')]
Ar_dist = distance.squareform(distance.pdist(df, metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
PCoA = skbio.stats.ordination.pcoa(DM_dist)
addedsamples = PCoA.samples.copy()
addedsamples.set_index(pd.MultiIndex.from_tuples(df.index, names=['Patient','time','type']),inplace=True)
addedsamples.reset_index(inplace=True)
addedsamples.time = addedsamples.time.astype('int')
oral_merged_metadata = pd.merge(merged_metadata, addedsamples, how='inner', left_on=['id','time_point'], right_on=['Patient', 'time'])
sns.scatterplot(data=oral_merged_metadata, x='PC1', y='PC2', hue='PPI', ax=oral_ax, palette='colorblind')

#plt.tight_layout()
#plt.savefig("results/PCoA.pdf")
plt.show()
