#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata for network annotation
Theo Portlock
'''
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import train_test_split

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_theme(font_scale=0.6)

variable = 'MELD'
taxaType='genus'
samples_metadata = pd.read_csv("newnewmetadata.csv").set_index('Sample ID').dropna(axis=0, subset=['MELD'])

gmsp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
gmsp_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
samples_gtaxonomy = gmsp_samples.join(gmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
gmspdf = samples_gtaxonomy.add_prefix('gut ')

omsp_samples = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)
omsp_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)
samples_otaxonomy = omsp_samples.join(omsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
omspdf = samples_otaxonomy.add_prefix('oral ')

df = omspdf.join(gmspdf, how='inner')
#df = df.replace(0,np.nan).dropna(axis=1, thresh=30)
#df = df.replace(np.nan,0)
df = df.join(samples_metadata[variable], how='inner').astype('float').dropna(axis=0)

plotdf = pd.DataFrame(index=df.drop(variable, axis=1).columns)
for i in range(1000):
    reg = ExtraTreesRegressor(n_estimators=100)
    X = df.drop(variable, axis=1)
    y = df.xs(variable, axis=1)
    X_train, X_test, y_train, y_test = train_test_split(X, y)
    reg.fit(X_train, y_train)
    plotdf[i] = reg.feature_importances_

cormat = df.corr(method='spearman').drop(variable)
cormat['imp'] = plotdf.mean(axis=1)
cormat = cormat[['imp', variable]]
cormat.Site=''
cormat.loc[cormat.index.str.contains('gut'), 'Site'] = 'gut'
cormat.loc[cormat.index.str.contains('oral'), 'Site'] = 'oral'
fig = sns.scatterplot(data=cormat, x=variable, y='imp', s=15)
sigcor = cormat[cormat.imp > cormat.imp.quantile(.95)][[variable, 'imp']]
[fig.annotate(species, values) for species, values in sigcor.iterrows()]
plt.xlim(-0.5,0.5)
fig.set(xlabel='MELD Spearman correlation', ylabel='Random forest feature importance')
plt.savefig('results/CorrelationvsImportanceMeld.pdf')
plt.show()

'''
# test individual
new = df.drop(variable, axis=1).join(samples_metadata, how='inner')
sns.swarmplot(data = new, x='Days after treatment', y='gut Atopobiaceae', hue='Type')
sns.boxplot(data = new, x='Days after treatment', y='gut Atopobiaceae', hue='Type')
'''

