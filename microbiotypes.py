#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Plotting an annotated seaborn correlation heatmap between MSP data and metadata
Theo Portlock
'''
%autoindent

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import skbio
from scipy.spatial import distance
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.cm as cm

msp_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)
#msp_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
samples_metadata = pd.read_csv('metadata.csv', index_col=0)
msp_samples = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)
#msp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)

# Join taxo information
taxaType = "species"
taxonomy_samples = msp_samples.join(msp_taxonomy[taxaType], how='inner').groupby(taxaType).sum()
samples_taxonomy = taxonomy_samples.T
samples_taxonomyMetadata = samples_taxonomy.join(samples_metadata, how='inner')

Ar_dist = distance.squareform(distance.pdist(samples_taxonomy, metric="braycurtis"))
#Ar_dist = distance.squareform(distance.jensenshannon(samples_taxonomy))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
PCoA = skbio.stats.ordination.pcoa(DM_dist)
addedsamples = PCoA.samples.copy()
addedsamples.set_index(samples_taxonomy.index,inplace=True)
gut_merged_metadata = addedsamples.join(samples_metadata)
#sns.scatterplot(data=gut_merged_metadata, x='PC1', y='PC2', hue='Aetiology', palette='colorblind')
#sns.scatterplot(data=gut_merged_metadata, x='PC1', y='PC2', hue='MELD', palette='coolwarm')

#distance.jensenshannon(a,b)
#plt.show()
#plt.tight_layout()
#plt.savefig("results/patric_metadata-clustermap.pdf")
'''
import plotly.express as px
merged_metadata = addedsamples.join(samples_metadata)

fig = px.scatter(
        merged_metadata.sort_values("Days after treatment"),
        x='PC1',
        y='PC2',
        color='Type',
        #size=80,
        #symbol='Aetiology',
        animation_frame='Days after treatment',
        #animation_group='Days after treatment',
        )
fig.show()

        

# elbow
Nc = range(1, 20)
kmeans = [KMeans(n_clusters=i) for i in Nc]
score = [kmeans[i].fit(addedsamples).score(addedsamples) for i in range(len(kmeans))]
plt.plot(Nc,score)
plt.xlabel('Number of Clusters')
plt.ylabel('Score')
plt.title('Elbow Curve')
plt.show()


X = addedsamples.iloc[:,:2].values

#range_n_clusters = [2, 3, 4, 5, 6]
range_n_clusters = [3]

for n_clusters in range_n_clusters:
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)
    ax1.set_xlim([-0.1, 1])
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)
    silhouette_avg = silhouette_score(X, cluster_labels)
    print("For n_clusters =", n_clusters,
          "The average silhouette_score is :", silhouette_avg)
    sample_silhouette_values = silhouette_samples(X, cluster_labels)
    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]
        ith_cluster_silhouette_values.sort()
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        y_lower = y_upper + 10  # 10 for the 0 samples
    ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")
    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
    # 2nd Plot showing the actual clusters formed
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                c=colors, edgecolor='k')
    # Labeling the clusters
    centers = clusterer.cluster_centers_
    # Draw white circles at cluster centers
    ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
                c="white", alpha=1, s=200, edgecolor='k')
    for i, c in enumerate(centers):
        ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
                    s=50, edgecolor='k')
    ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("Feature space for the 1st feature")
    ax2.set_ylabel("Feature space for the 2nd feature")
    plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                  "with n_clusters = %d" % n_clusters),
                 fontsize=14, fontweight='bold')
plt.show()

samples_metadata['cluster'] = cluster_labels
samples_taxonomyMetadata = samples_taxonomy.join(samples_metadata, how='inner')
samples_taxonomyMetadata['cluster'] = cluster_labels
sampmean = samples_taxonomyMetadata.groupby('cluster').mean()
fsm = sampmean.loc[:,samples_metadata.loc[:, samples_metadata.columns.isin(sampmean.columns)].columns]
nfsm = pd.DataFrame(scaler.fit_transform(fsm), columns=fsm.columns,index=fsm.index)
sns.clustermap(nfsm,cmap='coolwarm')

