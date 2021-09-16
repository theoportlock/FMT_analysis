#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import plotly.express as px
from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans

df = datasets.load_iris(as_frame=True).frame
#df = datasets.load_diabetes(as_frame=True).frame
scaledDf = StandardScaler().fit_transform(df)
pca = PCA()
results = pca.fit_transform(scaledDf)
df['PC1'], df['PC2'] = results[:,0], results[:,1]
n_clusters = range(2,16)
models = [KMeans(n_clusters=i).fit(df[['PC1', 'PC2']]) for i in n_clusters]
sscores = pd.Series([silhouette_score(df[['PC1', 'PC2']], i.labels_) for i in models], index=n_clusters)
print(sscores)
#df['Cluster'] = models[sscores.reset_index(drop=True).idxmax()].labels_
df['Cluster'] = models[1].labels_
sns.scatterplot(data = df, x='PC1', y='PC2', hue='Cluster', palette='colorblind')
stats = fdf.groupby('Days after treatment').apply(lambda group: mannwhitneyu(baseline[y], group[y])[1])

# need to look at all combinations of groups and do t test between them

df.groupby('Cluster').mean()
plt.show()
top=0.98,
bottom=0.045,
left=0.045,
right=0.955,
hspace=0.235,
wspace=0.17
