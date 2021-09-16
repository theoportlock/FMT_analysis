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
from matplotlib.patches import Ellipse

samples_metadata = pd.read_csv('newmergedmetadata.csv', index_col='ID_x')
samples_gutmsp = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0).T
samples_oralmsp = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0).T

variables = ['Days after treatment', 'Patient_Id', 'Type']
samples_mgutmsp = samples_gutmsp.join(samples_metadata[variables]).drop_duplicates()
samples_moralmsp = samples_oralmsp.join(samples_metadata[variables]).drop_duplicates()
samples_mgutmsp.loc[samples_mgutmsp.Type == 'DONOR', 'Days after treatment'] = 0

# Compute bray distances and plot
fig, (oral_ax, gut_ax) = plt.subplots(1, 2, sharey=True)
gut_ax.set_title("Gut")
oral_ax.set_title("Oral")

#df = StandardScaler().fit_transform(samples_mgutmsp.set_index(variables))
df = samples_gutmsp.values
Ar_dist = distance.squareform(distance.pdist(df, metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
PCoA = skbio.stats.ordination.pcoa(DM_dist, number_of_dimensions=2)
results = PCoA.samples.copy()
samples_mgutmsp['PC1'], samples_mgutmsp['PC2'] = results.iloc[:,0].values, results.iloc[:,1].values
samples_mgutmsp.reset_index().groupby("Patient_Id").plot(x="PC1", y="PC2", ax=gut_ax, legend=False, color='gray', linewidth=0.5, label='_nolegend_')
sns.scatterplot(data=samples_mgutmsp, x='PC1',y='PC2', style='Type', hue='Days after treatment', palette='coolwarm' ,ax=gut_ax)
#confidence_ellipse(samples_mgutmsp.loc[samples_mgutmsp.Type == 'FMT', 'PC1'], samples_mgutmsp.loc[samples_mgutmsp.Type == 'FMT', 'PC2'], gut_ax, n_std=1, label='FMT', alpha=0.5, facecolor='pink' ,edgecolor='purple', zorder=0)

#df = StandardScaler().fit_transform(samples_moralmsp.set_index(variables))
df = samples_oralmsp.values
Ar_dist = distance.squareform(distance.pdist(df, metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
PCoA = skbio.stats.ordination.pcoa(DM_dist, number_of_dimensions=2)
results = PCoA.samples.copy()
samples_moralmsp['PC1'], samples_moralmsp['PC2'] = results.iloc[:,0].values, results.iloc[:,1].values
samples_moralmsp.reset_index().groupby("Patient_Id").plot(x="PC1", y="PC2", ax=oral_ax, legend=False, color='gray', linewidth=0.5, label='_nolegend_')
sns.scatterplot(data=samples_moralmsp, x='PC1',y='PC2', style='Type', hue='Days after treatment', palette='coolwarm', ax=oral_ax)

sns.set(rc={'figure.figsize':(11.7,8.27)})
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.yaxis
oral_ax.legend([],[], frameon=False)
plt.tight_layout()
plt.savefig("results/newPCoA.pdf")
plt.show()
'''

from ellipsefunc import confidence_ellipse
fig, ax = plt.subplots()

mdf = samples_mgutmsp.groupby(['Type', 'Days after treatment']).mean()
sns.scatterplot(data=mdf, x='PC1',y='PC2', style='Days after treatment', hue='Type')

#mdf = samples_moralmsp.groupby(['Type', 'Days after treatment']).agg({'PC1':['mean', 'std'],'PC2':['mean', 'std']})
mdf = samples_moralmsp.groupby(['Type', 'Days after treatment']).mean()
smdf = samples_moralmsp.groupby(['Type', 'Days after treatment']).std()

sns.scatterplot(data=mdf, x='PC1',y='PC2', style='Days after treatment', hue='Type', ax=ax)
confidence_ellipse(mdf.PC1, mdf.PC2, ax, n_std=1, label='test', edgecolor='blue', linestyle=':')
'''
