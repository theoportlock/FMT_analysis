#!/usr/bin/env python
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

# Load
data = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv",index_col=0)
names = pd.read_table("../downstream_data/P15952_20_03_sample_info_stool.txt")
sample_type = pd.read_csv("../downstream_data/PROFITplaceboTab.csv", index_col=0)
msp_to_species = pd.read_csv("../downstream_data/hs_v2_sp.csv", index_col=0).set_index('msp_name').dropna()
card = pd.read_csv("../downstream_data/hs_10_4_igc2.CARD.tsv", sep='\t')

# Format
sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))
formatted_names = names['User ID'].str.split(expand=True)
formatted_names.loc[formatted_names[0]=="Donor", 2] = formatted_names.loc[formatted_names[0]=="Donor", 1]
formatted_names.loc[formatted_names[0]=="Donor", 0] = formatted_names.loc[formatted_names[0]=="Donor", 1]
formatted_names.iloc[:,2] = formatted_names.iloc[:,2].str.replace("Stool","0")
formatted_names['ID'] = names['NGI ID'].str.replace('P15952_','FMT')
formatted_names = pd.merge(formatted_names, sample_type, left_on=0, right_on="Patient", how='left')
formatted_names.fillna("DONOR", inplace=True)
card["genename"] = card.ORF_ID.str.split(expand=True)[0]
card["genename"] = card["genename"].str[0:-2]
data.columns = pd.MultiIndex.from_frame(formatted_names[[0,2,'Type','ID']], names=['Patient','Days after treatment', 'Type', 'ID'])
data.reset_index(inplace=True)
variable = 'gene_name'
repl = msp_to_species[variable].to_dict()
fdata = data.replace(repl)
fdata = fdata[(~fdata['index'].str.contains("msp_"))].reset_index()
mapped = fdata[fdata["index"].isin(card.genename)]

variable = 'Best_Hit_ARO'
#variable = 'AMR Gene Family'
#variable = 'Drug Class'
#variable = 'Resistance Mechanism'

repl = card.set_index('genename')[variable].to_dict()
mapped_variable = mapped.replace(repl)
mapped_variable = mapped_variable.drop("level_0",axis=1)
mapped_variable = mapped_variable.groupby("index").sum().T
result = mapped_variable.reset_index()
result['Days after treatment'] = result['Days after treatment'].astype('int')
fmtresult = result[result.Type != 'DONOR']

# Plot
fmtresult.loc[fmtresult.Type == 'FMT'].groupby(['Days after treatment']).mean().plot.bar(stacked=True)
#fmtresult.groupby(['Days after treatment', 'Type']).mean().plot.bar(stacked=True)
plt.ylabel("Normalized AMR gene abundance")
plt.tight_layout()
#plt.show()
plt.savefig('results/antibio.pdf')
