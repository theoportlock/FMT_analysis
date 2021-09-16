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

gutmsp_taxonomy = pd.read_csv("taxo.csv", index_col=0)
gutmsp_samples = pd.read_csv("msp.csv", index_col=0)
samples_metadata = pd.read_csv('metadata.csv', index_col=0)

# Join taxonomic information
taxaType = "species"
gutsamples_taxonomy = msp_samples.join(gutmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
scaled_samples_taxonomy = StandardScaler().fit_transform(samples_taxonomy)

# Compute PCoA
Ar_dist = distance.squareform(distance.pdist(scaled_samples_taxonomy, metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
PCoA = skbio.stats.ordination.pcoa(DM_dist, number_of_dimensions=2)
results = PCoA.samples.copy()
samples_taxonomy['PC1'], samples_taxonomy['PC2'] = results.iloc[:,0].values, results.iloc[:,1].values

# Plot metadata label (for you it would be before and after surgery)
samples_taxonomyMetadata = samples_taxonomy.join(samples_metadata, how='inner')
sns.scatterplot(data = samples_taxonomyMetadata, x='PC1', y='PC2', hue='Gender', palette='colorblind')
plt.show()

# Compute best cluster
n_clusters = range(2,20)
models = [KMeans(n_clusters=i).fit(samples_taxonomy[['PC1', 'PC2']]) for i in n_clusters]
sscores = pd.Series([silhouette_score(samples_taxonomy[['PC1', 'PC2']], i.labels_) for i in models], index=n_clusters)
print(sscores)
samples_taxonomy['Cluster'] = models[sscores.reset_index(drop=True).idxmax()].labels_

# Plot best cluster
sns.scatterplot(data = samples_taxonomy, x='PC1', y='PC2', hue='Cluster', palette='colorblind')
plt.show()





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
