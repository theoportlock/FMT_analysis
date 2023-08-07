#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import functions
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shap

meta = pd.read_csv("../data/meta.csv").set_index('Sample ID')
gutmsp = pd.read_csv("../data/gutmsp.csv", index_col=0).T
oralmsp = pd.read_csv("../data/oralmsp.csv", index_col=0).T
guttaxo = pd.read_csv("../data/guttaxo.csv", index_col=0)
oraltaxo = pd.read_csv("../data/oraltaxo.csv", index_col=0)

#msp=oralmsp.copy()
msp=gutmsp.copy()
taxo=guttaxo.copy()

#corr = msp.T.join(taxo["species"]).set_index("species").T.corr(method='spearman')
corr = msp.T.join(taxo["species"]).set_index("species").T.cov()
G = functions.network(corr, thresh = 0.7)
functions.clusterplot(G)
plt.tight_layout(); plt.show()
clust = functions.cluster(G)
functions.annotateplot(G, clust)

for i in clust.unique():
    ncor = corr.loc[clust.loc[clust == i].index, clust.loc[clust == i].index]
    G = functions.network(ncor, thresh = 0.5)
    functions.clusterplot(G)
    plt.show()

#functions.annotateplot(G, final.to_frame().join(taxo).sort_values(0).set_index("species")[0])
shaps = functions.rfr(msp.join(meta['MELDPLUS']).dropna(), 'MELDPLUS')

a = (
    shaps.to_frame().join(taxo["species"])
    .groupby("species")
    .sum()
)
b = pd.concat([clust, a], axis=1)

fig, (ax1, ax2) = plt.subplots(2,1)
sns.stripplot(data=b, x="Group", y=0, size=3, color="black", ax=ax1)
sns.boxplot(data=b, x="Group", y=0, showfliers=False, boxprops=dict(alpha=.5), ax=ax1)

vals = msp.T.join(taxo["species"]).groupby("species").sum().sum(axis=1)
c = pd.concat([clust, vals], axis=1).dropna()
ax2 = functions.abund()

df = c.pivot(columns="Group")[0].fillna(0)
unclass = df[df.index.str.contains("unclassified")].sum()
df = df[~df.index.str.contains("unclassified")]
df.loc['unclassified'] = unclass
if df.shape[0] > 20:
    df.loc['other'] = df.loc[df.T.sum().sort_values(ascending=False).iloc[21:].index].sum()
df = df.loc[df.T.sum().sort_values().tail(20).index]
norm = df.T.div(df.sum(axis=0), axis=0)
norm.plot(kind='bar',stacked=True, width=0.9, cmap='tab20', ylim=(0,1), ax=ax2)
ax2.legend(title='Taxonomy', bbox_to_anchor=(1.001, 1), loc='upper left', fontsize='small')
ax2.ylabel('Relative abundance')
ax2.setp(ax.get_xticklabels(), rotation=40, ha="right")

plt.tight_layout()
plt.show()
