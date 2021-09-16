#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import plotly.express as px
from sklearn import datasets
#from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
from itertools import combinations
from scipy.stats import mannwhitneyu

# Data location
samples_metadata = pd.read_csv('metadata.csv', index_col=0)
gmsp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
gmsp_gtaxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)

taxaType='species'
gtaxonomy_samples = gmsp_samples.join(gmsp_gtaxonomy[taxaType], how='inner').groupby(taxaType).sum()

var = 'Aetiology'
#df = gtaxonomy_samples.drop('FMT123', axis=1).T
df = gtaxonomy_samples.T.join(samples_metadata[var].dropna(), how='inner')
scaledDf = StandardScaler().fit_transform(df.drop(var, axis=1))
#pca = PCA()
#tmp = pca.fit(scaledDf, df[var])
results = PCA(n_components=2).fit(scaledDf, df[var]).transform(scaledDf)
df['PC1'], df['PC2'] = results[:,0], results[:,1]
n_clusters = range(2,16)
models = [KMeans(n_clusters=i).fit(df[['PC1', 'PC2']]) for i in n_clusters]
sscores = pd.Series([silhouette_score(df[['PC1', 'PC2']], i.labels_) for i in models], index=n_clusters)
print(sscores)
df['Cluster'] = models[sscores.reset_index(drop=True).idxmax()].labels_
#df['Cluster'] = models[7].labels_
sns.scatterplot(data = df, x='PC1', y='PC2', hue='Cluster', palette='colorblind')
#stats = df.groupby('Days after treatment').apply(lambda group: mannwhitneyu(baseline[y], group[y])[1])

# need to look at all combinations of groups and do t test between them
lut = dict(zip(df.Cluster.unique(), "rbg"))
row_colors = df.Cluster.map(lut)

addedcluster = samples_metadata.join(df['Cluster'])._get_numeric_data().dropna()
#addedcluster = df.join(df['Cluster'])

combs = list(combinations(addedcluster.Cluster.unique(), 2))
results = pd.DataFrame(index=addedcluster.columns)
for i, ii in combs:
    result = {}
    for j in addedcluster.columns:
        try:
            result[j] = mannwhitneyu(addedcluster.loc[addedcluster.Cluster == i,j].values, addedcluster.loc[addedcluster.Cluster == ii,j].values)[1]
        except:
            result[j] = np.nan 
    results[f'{i}_{ii}'] = result.values()

results < 0.05

sns.clustermap(df.xs(['PC1','PC2'], axis=1))
plt.show()

df2 = samples_metadata.join(df[['Cluster','PC1','PC2']])
fig, ax = plt.subplots()
sns.scatterplot(data = df2, x='PC1', y='PC2', style='Cluster', hue=var, size='Days after treatment', palette='colorblind', ax=ax)
#sns.kdeplot(data = df2.query("Type == 'FMT'"), x='PC1', y='PC2', cmap='Blues', ax=ax)
#sns.kdeplot(data = df2.query("Type == 'PLACEBO'"), x='PC1', y='PC2', shade=False, cmap='Reds', ax=ax)
#sns.scatterplot(data = df2, x='PC1', y='PC2', style='Cluster', hue='MELD', size='Days after treatment', palette='coolwarm')
sns.scatterplot(data = df2, x='PC1', y='PC2', style='Type', hue='MELD', size='Days after treatment', palette='coolwarm')

plt.show()
