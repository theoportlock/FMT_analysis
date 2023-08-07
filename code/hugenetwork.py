#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Plotting an annotated seaborn correlation heatmap between MSP data and metadata
Theo Portlock
'''

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

meta = pd.read_csv('../data/cleanMeta.csv').set_index('Sample ID').fillna(0)._get_numeric_data()
meta.columns = meta.columns.str.replace(' ','_')
meta = meta['MELDPLUS'].to_frame()
othermeta = pd.read_csv('../data/plasmaBileAcid.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('plasma ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/plasmaTryptophan.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('plasma ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/stoolBileAcid.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('stool ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/stoolNmr.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('stool ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/stoolTryptophan.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('stool ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/urineBileAcid.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('urine ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/urineNmr.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('urine ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/urineTryptophan.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('urine ')
meta = meta.join(othermeta)

meta = meta.loc[meta.MELDPLUS != 0]
meta.dropna(inplace=True)

taxo = pd.read_csv("../data/guttaxo.csv", index_col=0)
msp = pd.read_csv("../data/gutmsp.csv", index_col=0)
taxaType='species'
mspTaxo = msp.join(taxo[taxaType], how='inner').groupby(taxaType).sum().T
mspTaxo = mspTaxo.add_prefix('gut ')
meta = meta.join(mspTaxo, how='inner').astype('float')

taxo = pd.read_csv("../data/oraltaxo.csv", index_col=0)
msp = pd.read_csv("../data/oralmsp.csv", index_col=0)
taxaType='species'
mspTaxo = msp.join(taxo[taxaType], how='inner').groupby(taxaType).sum().T
mspTaxo = mspTaxo.add_prefix('oral ')
meta = meta.join(mspTaxo, how='inner').astype('float')

kegg = pd.read_csv("../data/gutkegg.csv", index_col=0).T
meta = meta.join(kegg).dropna()

patric = pd.read_csv("../data/gutpatric.csv", index_col=0).T
meta = meta.join(patric).dropna()

pfam = pd.read_csv("../data/gutpfam.csv", index_col=0).T
meta = meta.join(pfam).dropna()

card = pd.read_csv("../data/gutcard.csv", index_col=0).T
meta = meta.join(card).dropna()

cazyme = pd.read_csv("../data/gutcazy.csv", index_col=0).T
meta.drop(['CBM26', 'GH3', 'GT87', 'CBM27'], axis=1, inplace=True)
meta = meta.join(cazyme).dropna()

# Calculate and format correlations
correlationArray, uncorrectedPValueArray = spearmanr(meta)
correlations = pd.DataFrame(
    correlationArray,
    index=meta.columns,
    columns=meta.columns)
uncorrectedPValues = pd.DataFrame(
    uncorrectedPValueArray,
    index=meta.columns,
    columns=meta.columns)
filteredUncorrectedPValues = uncorrectedPValues.loc[
    (uncorrectedPValues.sum(axis=1) != 0),
    (uncorrectedPValues.sum(axis=0) != 0)]
filteredCorrelations = correlations.loc[
    (correlations.sum(axis=1) != 0),
    (correlations.sum(axis=0) != 0)]
significantMatrix = pd.DataFrame(
    multipletests(filteredUncorrectedPValues.values.flatten())[0].reshape(filteredUncorrectedPValues.shape),
    index = filteredUncorrectedPValues.index,
    columns = filteredUncorrectedPValues.columns)

cor = filteredCorrelations.stack().to_frame()
sig = significantMatrix.stack().to_frame()
fin = cor[sig].dropna()

uniq = pd.DataFrame(fin.reset_index().iloc[:, 0].unique())
uniq[0].apply(lambda x: ' '.join(x.split(' ')[1:]))
uniq['Site'] = uniq[0].str.split().str.get(0)
uniq['Rename'] = uniq[0].apply(lambda x: ' '.join(x.split(' ')[1:]))

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


var = 'MELDPLUS'
df = meta.copy().dropna()
df.loc[:] = StandardScaler().fit_transform(df)

i = df[var].unique()[0]
for i in df[var].unique():
    testdf = df.copy()
    model = ExtraTreesRegressor(n_estimators=500, n_jobs=-1,random_state=1)
    #model = ExtraTreesClassifier(n_estimators=500, n_jobs=-1,random_state=1)
    X = testdf.drop(var, axis=1)
    y = testdf.xs(var, axis=1)
    #y = pd.get_dummies(testdf.xs(var, axis=1))['DONOR'] * -1 + 1
    #X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 1, stratify=True)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 1)
    model.fit(X_train, y_train)
    #model = pickle.load(open('../../../../ATLAS/newatlas/randomforest/code/LC', 'rb'))
    y_pred = model.predict(X_test)
    #y_pred = model.predict(X)

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

0   K07000
1   K11267
2   K03925
3   K07169
4   K13929
5   K10560
6   K11136
7   K15449
8   K17881
9   K19183
10  K19349
11  K00942
12  K18457
13  K20307
14  K14481
15  K08984
16  K09690
17  K02463
18  K03247
19  K12818
20  K16296
21  K17105
22  K13798
23  K01128
24  K01193
25  K19174
26  K07235
27  K11835
28  K14696
29  K18608

stool ba

                                               index         0
0                         Taurochenodeoxycholic Acid  0.012599
1                                Hyodeoxycholic Acid  0.012860
2  8(14),(5-beta)-Cholenic Acid-3-alpha, 12-alpha...  0.013461
3                             Glycoursocholanic Acid  0.014127
4   5-beta-Cholanic Acid-3-alpha, 6-alpha-diol-7-one  0.018781
5                          Glycoursodeoxycholic Acid  0.022287
6                          Tauroursodeoxycholic Acid  0.024634
7                                   Taurocholic Acid  0.025794
8                    Chenodeoxycholic Acid-3-Sulfate  0.074775
9                              Cholic Acid-3-Sulfate  0.075468
