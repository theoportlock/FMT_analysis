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

variable = 'MELD'
taxaType='species'
samples_metadata = pd.read_csv('metadata.csv', index_col=0).dropna(axis=0)

gmsp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
gmsp_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
samples_gtaxonomy = gmsp_samples.join(gmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
gmspdf = samples_gtaxonomy.add_prefix('gut ')

omsp_samples = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)
omsp_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)
samples_otaxonomy = omsp_samples.join(omsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
omspdf = samples_otaxonomy.add_prefix('oral ')

df = omspdf.join(gmspdf, how='inner')
df = df.join(samples_metadata[variable], how='inner').astype('float')

reg = ExtraTreesRegressor(n_estimators=10000)
X = df.drop(variable, axis=1)
y = df.xs(variable, axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, y)
reg.fit(X_train, y_train)
cormat = df.corr(method='spearman').drop('MELD')
cormat['imp'] = reg.feature_importances_
fig = px.scatter(cormat.reset_index(), x='MELD', y='imp', hover_name='index')
fig.show()
cormat[cormat.imp > cormat.imp.quantile(.95)][['MELD', 'imp']]
#cormat[['MELD', 'imp']].to_csv('randforimp.csv')
