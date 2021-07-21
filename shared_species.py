#!/usr/bin/env python
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

# Data location
samples_metadata = pd.read_csv('metadata.csv', index_col=0)
gmsp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
omsp_samples = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)
gmsp_gtaxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
omsp_otaxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)

taxaType='species'
gtaxonomy_samples = gmsp_samples.join(gmsp_gtaxonomy[taxaType], how='inner').groupby(taxaType).sum()
samples_gtaxonomy = gtaxonomy_samples.T
samples_taxonomyMetadata = samples_gtaxonomy.join(samples_metadata, how='inner')

samples_gutsum['Site'] = 'Gut'
samples_oralsum['Site'] = 'Oral'
joined = pd.concat([samples_gutsum, samples_oralsum]).reset_index()

#stopped here

# Format data

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
