#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
from scipy import stats
from scipy.spatial import distance
from scipy.stats import spearmanr
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import (accuracy_score, confusion_matrix, classification_report)
from sklearn.metrics import r2_score
from sklearn.metrics import explained_variance_score
from sklearn.metrics import plot_confusion_matrix
from sklearn.metrics import plot_precision_recall_curve
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shap
import skbio

#taxaType = 'genus'
#variable = 'MELD'
#var = variable
samples_metadata = pd.read_csv('../data/newnewmetadata.csv').set_index('Sample ID') 
#samples_metadata = samples_metadata.dropna(subset=variable)

msp_samples = pd.read_csv("../data/gutmsp.csv", index_col=0)
msp_taxonomy = pd.read_csv("../data/guttaxo.csv", index_col=0)
samples_taxonomy = msp_samples.reindex(msp_taxonomy.index, fill_value=0).T
taxonomy_samples = msp_samples.join(msp_taxonomy[taxaType], how='inner').groupby(taxaType).sum()
samples_taxonomy = taxonomy_samples.T
samples_taxonomy.loc[:] = StandardScaler().fit_transform(samples_taxonomy)
#mspdf = samples_taxonomy.join(samples_metadata[variable], how='inner').astype('float')
mspdf = samples_taxonomy.join(samples_metadata['Type'], how='inner')
#mspdf = samples_taxonomy.join(samples_metadata[variable], how='inner')

oralmsp_samples = pd.read_csv("../data/oralmsp.csv", index_col=0)
oralmsp_taxonomy = pd.read_csv("../data/oraltaxo.csv", index_col=0)
samples_oraltaxonomy = oralmsp_samples.join(oralmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
samples_oraltaxonomy.loc[:] = StandardScaler().fit_transform(samples_oraltaxonomy)
omspdf = samples_oraltaxonomy.join(samples_metadata[variable], how='inner').astype('float')

samples_kegg = pd.read_csv("../data/gutkegg.csv", index_col=0).T
samples_kegg.loc[:] = StandardScaler().fit_transform(samples_kegg)
kmdf = samples_kegg.join(samples_metadata[variable]).dropna()

samples_patric = pd.read_csv("gutpatric.csv", index_col=0).T
samples_patric.loc[:] = StandardScaler().fit_transform(samples_patric)
padf = samples_patric.join(samples_metadata[variable]).dropna()

samples_pfam = pd.read_csv("../data/gutpfam.csv", index_col=0).T
samples_pfam.loc[:] = StandardScaler().fit_transform(samples_pfam)
pfdf = samples_pfam.join(samples_metadata[variable]).dropna()

samples_card = pd.read_csv("../data/gutcard.csv", index_col=0).T
samples_card.loc[:] = StandardScaler().fit_transform(samples_card)
cadf = samples_card.join(samples_metadata[variable]).dropna()

samples_cazyme = pd.read_csv("../data/gutcazy.csv", index_col=0).T
samples_cazyme.loc[:] = StandardScaler().fit_transform(samples_cazyme)
czdf = samples_cazyme.join(samples_metadata[variable]).dropna()

samples_metab = pd.read_csv("../data/formatmetabolite.csv").dropna().set_index('Sample ID')
samples_metab.loc[:] = StandardScaler().fit_transform(samples_metab)
medf = samples_metab.join(samples_metadata[variable]).dropna()

samples_meta = samples_metadata._get_numeric_data()
samples_meta.loc[:] = StandardScaler().fit_transform(samples_meta)
medf = samples_meta.dropna(axis=1)


feature_imp = pd.DataFrame()
scorecurve = pd.DataFrame(columns=['scores', 'curves'])
scores = pd.Series()
curves = pd.Series()
shaps = pd.DataFrame()
plt.rcParams["figure.figsize"] = (7,7)

i = df[var].unique()[0]
for i in df[var].unique():
    #testdf = medf.copy()
    testdf = mspdf.copy()
    #model = ExtraTreesRegressor(n_estimators=500, n_jobs=-1,random_state=1)
    model = ExtraTreesClassifier(n_estimators=500, n_jobs=-1,random_state=1)
    X = testdf.drop(var, axis=1)
    #y = testdf.xs(var, axis=1)
    y = pd.get_dummies(testdf.xs(var, axis=1))['DONOR'] * -1 + 1
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 1, stratify=True)
    model.fit(X_train, y_train)
    model = pickle.load(open('../../../../ATLAS/newatlas/randomforest/code/LC', 'rb'))
    #y_pred = model.predict(X_test)
    y_pred = model.predict(X)
plot_roc_curve(model, X, y, pos_label=1)
plt.show()
plot_confusion_matrix(model, X, y, display_labels=['H', 'D'], colorbar=False, cmap='Reds')
plt.show()
plot_precision_recall_curve(model, X, y, pos_label=1)
plt.show()


    r2_score(y_true=y_test, y_pred=y_pred)
    explained_variance_score(y_true=y_test, y_pred=y_pred)

    feature_imp[i] = pd.Series(model.feature_importances_,index=X.columns)
    explainer = shap.TreeExplainer(model)
    #shaps[i] = pd.Series(explainer(X).values.sum(axis=0)[:,0], index=X.columns)
    shaps_values = explainer(X)
    meanabsshap = pd.Series( np.abs(shaps_values.values).mean(axis=0), index=X.columns)
    corrs = [spearmanr(shaps_values.values[:,x], X.iloc[:,x])[0] for x in range(len(X.columns))]
    final = meanabsshap * np.sign(corrs)
    final.fillna(0, inplace=True)
    shap.summary_plot(shaps_values, X)
    #shaps[i] = final

'''
sscores = scores.apply(lambda x: "{:.4f}".format(x))
shaps.columns = scores.index.str.cat(sscores, sep=" ")
'''
#plotdf = feature_imp[(stats.zscore(feature_imp) > 4).any(axis=1)]
#plotdf = feature_imp[(feature_imp > 0.02).any(axis=1)]
#sns.clustermap(plotdf, yticklabels=True, cmap='Reds')
#plotdf = shaps[(np.abs(stats.zscore(shaps)) > 6).any(axis=1)]
plotdf = shaps[(np.abs(shaps) > 0.0125).any(axis=1)]
plotdf = plotdf[~plotdf.index.str.contains('unclassified')]
#plotdf = shaps[(np.abs(stats.zscore(shaps, axis=1)) > 4.5).any(axis=1)]
#plotdf = shaps[(np.abs(shaps) > 0.5).any(axis=1)]
#sns.clustermap(plotdf, yticklabels=True, cmap='coolwarm', center=0, vmax=1, vmin=-1)
#sns.clustermap(plotdf, yticklabels=True, cmap='coolwarm', center=0, z_score=True,square=True)
#sns.heatmap(plotdf.apply(stats.zscore), yticklabels=True, cmap='coolwarm', center=0, square=True, xticklabels=True)
#sns.heatmap(plotdf, yticklabels=True, cmap='coolwarm', center=0, square=True, xticklabels=True)
#plt.savefig('../results/shap2.svg')
#plt.show()
#plotdf.join(taxo.set_index('species')['gp']).set_index('gp').to_csv('shap2.txt',sep=' ')

enrich = pd.read_csv('../data/S3enrich.tsv', sep= '\t',index_col=0)
enrich.loc[enrich.Enriched == 'Depleted', enrich.select_dtypes(include=['number']).columns] *= -1
enrich.columns = enrich.columns.str.replace(':.*', '', regex=True)
#enrich = enrich.reset_index().set_index(['msp', 'Unnamed', 'Enriched'])
#enrich.xs(enrich.Enriched == 'Enriched', level=0)
enrich = enrich.set_index('Unnamed').drop('Enriched', axis=1)
enrich = enrich[~enrich.index.str.contains('unclassified')]
enrich = enrich.T.groupby(enrich.columns).mean().T
enrich = enrich.groupby(enrich.index).mean()
enrich.rename({'ME/CFS': 'ME_CFS'}, axis=1, inplace=True)
#filt = enrich.loc[:,enrich.columns.isin(plotdf.columns)]
filt = enrich.loc[enrich.index.isin(plotdf.index),enrich.columns.isin(scores.index)]

#eplotdf = enrich[(np.abs(stats.zscore(enrich)) > 4.7).any(axis=1)]
#sns.heatmap(eplotdf, xticklabels=True, yticklabels=True, cmap='vlag', center=0, square=True)
plotdf = plotdf.reindex(sorted(plotdf.columns), axis=1)
filt = filt.reindex(sorted(filt.columns), axis=1)

#sns.clustermap(filt, yticklabels=True, cmap='coolwarm', center=0, z_score=True, row_cluster=False, col_cluster=False, vmax=4)
#sns.heatmap(filt.apply(stats.zscore).drop('ob', axis=1), yticklabels=True, cmap='coolwarm', center=0, square=True, xticklabels=True)
sns.heatmap(filt, yticklabels=True, cmap='coolwarm', center=0, square=True, xticklabels=True)
plt.savefig('../results/filt2.svg')
#sns.clustermap((plotdf*filt).dropna(axis=1), yticklabels=True, cmap='coolwarm', center=0, z_score=True)

sns.heatmap(plotdf, yticklabels=True, cmap='coolwarm', center=0, square=True, xticklabels=True)
ax = sns.clustermap(plotdf, yticklabels=True, cmap='coolwarm', center=0,xticklabels=True)
plt.savefig('../results/shap2.svg')
plt.show()


'''
K03852
K03690
K10041
K06877
K04652
K13369
K03436
K00889
K02487
K06887
K00619
K06610
K16168
K07491
K02855
K06175
K04733
K07701
K02026
K12058
'''

'''
SUCCESSFULL PREDICTIONS 
P0004
P0015
P0016
P0020
P0024
