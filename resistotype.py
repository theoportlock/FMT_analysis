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

samples_metadata = pd.read_csv('newmergedmetadata.csv').drop_duplicates(subset='ID_x', keep='first').set_index('ID_x').sort_index().dropna(thresh=20, axis=1)

card_samples = pd.read_csv("fullcardNorm.csv", index_col=0)
var = 'AMR Gene Family'
samples_card = card_samples.groupby(var).sum().T

variables = ['Days after treatment', 'Patient_Id', 'Type']
samples_mgutmsp = samples_metadata[variables].join(samples_card, how='inner')
samples_mgutmsp.loc[samples_mgutmsp.Type == 'DONOR', 'Days after treatment'] = 0

fig, ax = plt.subplots()
df = samples_mgutmsp._get_numeric_data().values
Ar_dist = distance.squareform(distance.pdist(df, metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
PCoA = skbio.stats.ordination.pcoa(DM_dist, number_of_dimensions=2)
results = PCoA.samples.copy()
samples_mgutmsp['PC1'], samples_mgutmsp['PC2'] = results.iloc[:,0].values, results.iloc[:,1].values
samples_mgutmsp.reset_index().groupby("Patient_Id").plot(x="PC1", y="PC2",legend=False, color='gray', linewidth=0.5, label='_nolegend_', ax=ax)
sns.scatterplot(data=samples_mgutmsp, x='PC1',y='PC2', style='Type', hue='Days after treatment', palette='coolwarm', ax=ax) 
#confidence_ellipse(samples_mgutmsp.loc[samples_mgutmsp.Type == 'FMT', 'PC1'], samples_mgutmsp.loc[samples_mgutmsp.Type == 'FMT', 'PC2'], gut_ax, n_std=1, label='FMT', alpha=0.5, facecolor='pink' ,edgecolor='purple', zorder=0)

sns.set(rc={'figure.figsize':(11.7,8.27)})
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.tight_layout()
plt.savefig("results/resistotype.pdf")
#plt.show()
