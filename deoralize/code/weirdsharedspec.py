#!/usr/bin/env python
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

samples_metadata = pd.read_csv('../../data/newnewmetadata.csv').dropna(subset=['Sample ID']).set_index('Sample ID').dropna(thresh=20, axis=1)
samples_metadata.loc[samples_metadata.Type == 'DONOR', 'Days after treatment'] = 0

gmsp_samples = pd.read_csv("../../data/gutmsp.csv", index_col=0)
omsp_samples = pd.read_csv("../../data/oralmsp.csv", index_col=0)
gmsp_gtaxonomy = pd.read_csv("../../data/guttaxo.csv", index_col=0)
omsp_otaxonomy = pd.read_csv("../../data/oraltaxo.csv", index_col=0)

taxaType='species'
v = ['Days after treatment','Type','id']
gtaxonomy_samples = gmsp_samples.join(gmsp_gtaxonomy[taxaType], how='inner').groupby(taxaType).sum() > 0
samples_gtaxonomy = gtaxonomy_samples.T

otaxonomy_samples = omsp_samples.join(omsp_otaxonomy[taxaType], how='inner').groupby(taxaType).sum() > 0
samples_otaxonomy = otaxonomy_samples.T

sharedind = pd.Index.intersection(samples_gtaxonomy.index, samples_otaxonomy.index)
sharedcol = pd.Index.intersection(samples_gtaxonomy.columns, samples_otaxonomy.columns)

sharedvals = samples_gtaxonomy.loc[sharedind, sharedcol] & samples_otaxonomy.loc[sharedind, sharedcol]
sharedvals_metadata = samples_otaxonomy.join(samples_metadata[['Days after treatment','Type','id']], how='inner')
gvals = sharedvals_metadata.reset_index().groupby('id').count()['Days after treatment']
gfvals = gvals.loc[(gvals == 4) | (gvals == 1)].index
df = sharedvals_metadata[sharedvals_metadata['id'].isin(gfvals)].set_index(v)

df = df.droplevel(2).groupby(level=[0,1]).mean().ge(0.8)
df = df.sort_index(level=1)
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set(font_scale=0.6)
plotdf = df.T[df.sum() != 0]
plotdf = plotdf.loc[~plotdf.index.str.contains('unclassified')]
#plotdf = plotdf.loc[plotdf.xs(0, level=0, axis=1).any(axis=1)]
g = sns.heatmap(plotdf, yticklabels=True, cbar=False)
g.set_xticklabels(g.get_xticklabels(), rotation = 90)
plt.tight_layout()
#plt.savefig('results/sharedspecies.pdf')
plt.show()
'''
#test
individ = sharedvals_metadata.set_index(['Days after treatment', 'Type', 'id']).sum(axis=
1)
individ
individ.to_frame()
individ.to_frame().reset_index()
inddf = individ.to_frame().reset_index()
inddf
sns.boxplot(data=inddf, x='Days after treatment', y=0, hue='id', style='Type')
sns.stripplot(data=inddf, x='Days after treatment', y=0, hue='id', style='Type')
sns.stripplot(data=inddf, x='Days after treatment', y=0, hue='id', edgecolor='Type')
sns.stripplot(data=inddf, x='Days after treatment', y=0, hue='id')
plt.show()
sns.stripplot(data=inddf, x='Days after treatment', y=0, hue='id')
plt.show()
sns.pointplot(data=inddf, x='Days after treatment', y=0, hue='id')
plt.show()
history
'''
