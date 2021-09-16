#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
%autoindent
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor

variable = 'MELD'
samples_metadata = pd.read_csv('metadata.csv', index_col=0).dropna(axis=0)

msp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
#msp_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
#taxaType='species'
#taxonomy_samples = msp_samples.join(msp_taxonomy[taxaType], how='inner').groupby(taxaType).sum()
#samples_taxonomy = taxonomy_samples.T
oralmsp_samples = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)
#oralmsp_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)
#samples_oraltaxonomy = oralmsp_samples.join(oralmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T

omspdf = samples_oraltaxonomy.join(samples_metadata[variable], how='inner').astype('float')
mspdf = samples_taxonomy.join(samples_metadata[variable], how='inner').astype('float')
smdf = pd.read_csv("antismashNorm.csv", index_col=0).T.join(samples_metadata[variable]).dropna()
padf = pd.read_csv("patricNorm.csv", index_col=0).T.join(samples_metadata[variable]).dropna()
pfdf = pd.read_csv("pfamNorm.csv", index_col=0).T.join(samples_metadata[variable]).dropna()
cadf = pd.read_csv("cardNorm.csv", index_col=0).T.join(samples_metadata[variable]).dropna()
czdf = pd.read_csv("cazymeNorm.csv", index_col=0).T.join(samples_metadata[variable]).dropna()

def runmodel(df):
    #reg = RandomForestRegressor()
    cv = LeaveOneOut()
    model = ExtraTreesRegressor(random_state=1)
    X = df.drop(variable, axis=1)
    y = df.xs(variable, axis=1)
    scores = cross_val_score(model, X, y, scoring='neg_mean_absolute_error', cv=cv, n_jobs=-1)
    scores = absolute(scores)
    return scores

resultdf = pd.DataFrame()
resultdf['Gut Species'] = runmodel(mspdf)
resultdf['Oral Species'] = runmodel(omspdf)
resultdf['Secondary Metabolites'] = runmodel(smdf)
resultdf['Virulence Factors'] = runmodel(padf)
resultdf['Protein Family'] = runmodel(pfdf)
resultdf['Antimircobial Resistance'] = runmodel(cadf)
resultdf['Cazymes'] = runmodel(czdf)

plotdf.T.melt()
sns.pointplot(data=resultdf.T.melt(),x='variable', y='value')
plt.xlabel('Model Accuracy (%)')
plt.xlim(0,100)
plt.ylabel('Dataset used for MELD prediction')
plt.show()
plt.tight_layout()
#plt.savefig('results/OvG_randomfores_MELD.pdf')
