#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
%autoindent
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial import distance
import skbio

def evaluate(model, test_features, test_labels):
    predictions = model.predict(test_features)
    errors = abs(predictions - test_labels)
    mape = 100 * np.mean(errors / test_labels)
    accuracy = 100 - mape
    print('Model Performance')
    print('Average Error: {:0.4f} degrees.'.format(np.mean(errors)))
    print('Accuracy = {:0.2f}%.'.format(accuracy))
    return accuracy

variable = 'MELD'
samples_metadata = pd.read_csv('metadata.csv', index_col=0).dropna(axis=0)
msp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
msp_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
taxaType='species'
taxonomy_samples = msp_samples.join(msp_taxonomy[taxaType], how='inner').groupby(taxaType).sum()
samples_taxonomy = taxonomy_samples.T

oralmsp_samples = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)
oralmsp_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)
samples_oraltaxonomy = oralmsp_samples.join(oralmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
omspdf = samples_oraltaxonomy.join(samples_metadata[variable], how='inner').astype('float')

mspdf = samples_taxonomy.join(samples_metadata[variable], how='inner').astype('float')
#df = pd.get_dummies(samples_metadata)

smdf = pd.read_csv("antismashNorm.csv", index_col=0).T.join(samples_metadata[variable]).dropna()
padf = pd.read_csv("patricNorm.csv", index_col=0).T.join(samples_metadata[variable]).dropna()
pfdf = pd.read_csv("pfamNorm.csv", index_col=0).T.join(samples_metadata[variable]).dropna()
cadf = pd.read_csv("cardNorm.csv", index_col=0).T.join(samples_metadata[variable]).dropna()
czdf = pd.read_csv("cazymeNorm.csv", index_col=0).T.join(samples_metadata[variable]).dropna()

smdf = pd.read_csv("antismashNorm.csv", index_col=0).T
padf = pd.read_csv("patricNorm.csv", index_col=0).T
pfdf = pd.read_csv("pfamNorm.csv", index_col=0).T
cadf = pd.read_csv("cardNorm.csv", index_col=0).T
czdf = pd.read_csv("cazymeNorm.csv", index_col=0).T

def runmodel(df):
    #reg = RandomForestRegressor()
    reg = ExtraTreesRegressor()
    X = df.drop(variable, axis=1)
    y = df.xs(variable, axis=1)
    X_train, X_test, y_train, y_test = train_test_split(X, y)
    reg.fit(X_train, y_train)
    return evaluate(reg, X_test, y_test)

plotdf = pd.DataFrame()

for i in range(10):
    final = pd.Series()
    final['Gut Species'] = runmodel(mspdf)
    final['Oral Species'] = runmodel(omspdf)
    final['Secondary Metabolites'] = runmodel(smdf)
    final['Virulence Factors'] = runmodel(padf)
    final['Protein Family'] = runmodel(pfdf)
    final['Antimircobial Resistance'] = runmodel(cadf)
    final['Cazymes'] = runmodel(czdf)
    final.sort_values(ascending=False, inplace=True)
    plotdf[i] = final

plotdf.T.melt()
sns.pointplot(data=plotdf.T.melt(),x='variable', y='value')
plt.xlabel('Model Accuracy (%)')
plt.xlim(0,100)
plt.ylabel('Dataset used for MELD prediction')
plt.show()
plt.tight_layout()
#plt.savefig('results/OvG_randomfores_MELD.pdf')

df = mspdf.copy()
reg = ExtraTreesRegressor()
X = df.drop(variable, axis=1)
y = df.xs(variable, axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, y)
reg.fit(X_train, y_train)
cormat = mspdf.corr(method='spearman').drop('MELD')
cormat['imp'] = reg.feature_importances_
fig = px.scatter(cormat.reset_index(), x='MELD', y='imp', hover_name='index')
fig.show()
