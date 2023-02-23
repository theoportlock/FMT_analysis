#!/usr/bin/env python
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy.stats import mannwhitneyu

samples_metadata = pd.read_csv('newnewmetadata.csv').dropna(subset=['Sample ID']).set_index('Sample ID').dropna(thresh=20, axis=1)
samples_metadata.loc[samples_metadata.Type == 'DONOR', 'Days after treatment'] = 0

gmsp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
omsp_samples = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)
gmsp_gtaxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
omsp_otaxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)

taxaType='genus'
gtaxonomy_samples = gmsp_samples.join(gmsp_gtaxonomy[taxaType], how='inner').groupby(taxaType).sum() > 0
samples_gtaxonomy = gtaxonomy_samples.T

variables=['id','Donor', 'Type', 'Days after treatment']
df = samples_gtaxonomy.join(samples_metadata[variables], how='inner')
donordf=df.loc[df.Type == 'DONOR']
patientdf=df.loc[(df.Type != 'DONOR') & (df.Type == 'FMT')]
        
#shareddf = patientdf.drop(variables,axis=1) & donordf.drop(variables,axis=1).loc[donordf.id == patientdf.Donor]
shareddf = pd.DataFrame(index=patientdf.drop(variables, axis=1).columns)
for i in patientdf.index:
    #shareddf[i] = (patientdf.loc[i].drop(variables) & donordf.loc[donordf.id == patientdf.loc[i,'Donor']])
    shareddf[i] = (donordf.loc[donordf.id == patientdf.loc[i,'Donor']].drop(variables,axis=1) & patientdf.loc[i].drop(variables)).squeeze()

sharedvals_metadata = shareddf.T.join(samples_metadata[variables])

#otaxonomy_samples = omsp_samples.join(omsp_otaxonomy[taxaType], how='inner').groupby(taxaType).sum() > 0
#samples_otaxonomy = otaxonomy_samples.T

#sharedind = pd.Index.intersection(samples_gtaxonomy.index, samples_otaxonomy.index)
#sharedcol = pd.Index.intersection(samples_gtaxonomy.columns, samples_otaxonomy.columns)

#sharedvals = samples_gtaxonomy.loc[sharedind, sharedcol] & samples_otaxonomy.loc[sharedind, sharedcol]
#sharedvals_metadata = samples_otaxonomy.join(samples_metadata[['Days after treatment','Type','id']], how='inner')
#sharedvals_metadata = samples_otaxonomy.join(samples_metadata[['Days after treatment','Type','id','Donor']], how='inner')

df = sharedvals_metadata.drop('id',axis=1).groupby(['Days after treatment','Type']).mean().ge(0.9)
df = df.sort_index(level=1)
sns.set(rc={'figure.figsize':(14.0,8.27)})
sns.set(font_scale=0.6)
g = sns.heatmap(df.T[df.sum() != 0], yticklabels=True, cbar=False)
g.set_xticklabels(g.get_xticklabels(), rotation = 0)
plt.tight_layout()
plt.savefig('results/engraphment.pdf')
plt.show()
