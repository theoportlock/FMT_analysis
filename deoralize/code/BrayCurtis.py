#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Plotting an MSP PCoA 
'''

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import skbio
from scipy.spatial import distance
from scipy.stats import mannwhitneyu
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

taxo = pd.read_csv('../../data/guttaxo.csv', index_col=0)
msp = pd.read_csv('../../data/gutmsp.csv', index_col=0)
meta = pd.read_csv('../../data/cleanMeta.csv', index_col=0)

# Join taxonomic information
taxaType = "species"
mspTaxo = msp.join(taxo[taxaType], how='inner').groupby(taxaType).sum().T
scaledMspTaxo = StandardScaler().fit_transform(mspTaxo)

# Compute PCoA
Ar_dist = distance.squareform(distance.pdist(scaledMspTaxo, metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
PCoA = skbio.stats.ordination.pcoa(DM_dist, number_of_dimensions=2)
results = PCoA.samples.copy()
mspTaxo['PC1'], mspTaxo['PC2'] = results.iloc[:,0].values, results.iloc[:,1].values

# Plot metadata label
mspTaxoMeta = mspTaxo.join(meta, how='inner')
sns.scatterplot(data=mspTaxoMeta, x='PC1', y='PC2', hue='Gender', palette='colorblind')
plt.show()

# Compute best cluster
n_clusters = range(2,20)
models = [KMeans(n_clusters=i).fit(mspTaxoMeta[['PC1', 'PC2']]) for i in n_clusters]
sscores = pd.Series([silhouette_score(mspTaxoMeta[['PC1', 'PC2']], i.labels_) for i in models], index=n_clusters)
print(sscores)
mspTaxoMeta['Cluster'] = models[sscores.reset_index(drop=True).idxmax()].labels_

# Plot best cluster
sns.scatterplot(data = mspTaxoMeta, x='PC1', y='PC2', hue='Cluster', palette='colorblind')
plt.show()
