#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Plotting an MSP PCoA 
'''

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import skbio
from scipy.spatial import distance
from scipy.stats import mannwhitneyu
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

meta = pd.read_csv("../../data/newnewmetadata.csv").set_index('Sample ID').fillna(0)
msp = pd.read_csv("../../data/gutmsp.csv", index_col=0)
taxaType='species'
pval=0.05
taxo = pd.read_csv("../../data/guttaxo.csv", index_col=0)
msptaxo = msp.join(taxo[taxaType], how='inner').groupby(taxaType).sum().T
msptaxo = msptaxo[msptaxo.columns[~msptaxo.columns.str.contains('unclassified')]]
msptaxo.columns = msptaxo.columns.str.replace(r' /.*', '', regex=True)
msptaxo.columns = msptaxo.columns.str.replace(r' &.*', '', regex=True)
#msptaxo = msptaxo.loc[:, msptaxo.nunique() > 20]
v=['Days after treatment', 'Type', 'id', 'Donor']
df = msptaxo.join(meta[v], how='inner').set_index(v)
scaled_samples_taxonomy = StandardScaler().fit_transform(df)

Ar_dist = distance.squareform(distance.pdist(scaled_samples_taxonomy, metric="braycurtis"))
dist = pd.DataFrame(Ar_dist, index=df.index, columns=df.index)

sns.clustermap(data = dist,xticklabels=True, yticklabels=True, cmap='binary')

donordf=dist.xs('DONOR', level=1).droplevel([0,2])
patientdf=dist.drop('DONOR',level=1).xs('FMT', level=1)
a = patientdf.xs('DONOR', level=1, axis=1).droplevel([0,2], axis=1)
        
b = a.reset_index()
b['Beta diversity'] = 0
for i in b.index:
    b.loc[i, 'Beta diversity'] = b.loc[i, b.loc[i, 'Donor']]

fig, ax = plt.subplots(1,1)
plt.rcParams["figure.figsize"] = (5,5)
sns.pointplot(data =b, x='Days after treatment', y='Beta diversity', hue='id', ax=ax, scale=0.5)
sns.boxplot(data=b, y='Beta diversity', x='Days after treatment', linewidth=1, ax=ax)
#sns.stripplot(data=b, y='Beta diversity', x='Days after treatment', linewidth=1, ax=ax)

df = b
x='Days after treatment'
y='Beta diversity'
import itertools
# Stats
box_pairs = list(itertools.combinations(df[x].unique(),2))
test_short_name = 'M.W.W'
pvalues = {}
for pair in box_pairs:
    data1 = df.groupby([x])[y].get_group(pair[0])
    data2 = df.groupby([x])[y].get_group(pair[1])
    stat, p = mannwhitneyu(data1, data2)
    pvalues[pair] = p
pvals = pd.Series(pvalues)
qvals = pd.Series(fdrcorrection(pvals)[1], index=pvals.index)
sigqvals = qvals.loc[qvals < 0.05]

if len(sigqvals) > 0:
    ax, test_results = add_stat_annotation(
        ax,
        box_pairs=sigqvals.index,
        data=df,
        x=x,
        y=y,
        perform_stat_test=False,
        pvalues=sigqvals,
        test_short_name='M.W.W',
        text_format='full',
        #loc='outside',
        #verbose=2,
        )

plt.legend(title=y, bbox_to_anchor=(1.001, 1), loc='upper left', fontsize='small')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'../results/{y}.svg')
plt.show()

'''
# betadiv doesnt corr very well with anything including meld
engraphment directly doesnt correlate with meld, oral disruption does!
mergedagain = b[['Days after treatment', 'id', 'Beta diversity']].merge(meta, on=['Days after treatment', 'id'], how='inner')
#sns.clustermap(mergedagain._get_numeric_data().corr().fillna(0),xticklabels=True, yticklabels=True)
fig, ax = plt.subplots(1,1)
plt.rcParams["figure.figsize"] = (5,5)
mergedagain = mergedagain.loc[mergedagain.MELD != 0]
sns.pointplot(data =mergedagain, x='Days after treatment', y='MELD', hue='id', ax=ax, scale=0.5)
sns.boxplot(data=mergedagain, y='MELD', x='Days after treatment', linewidth=1, ax=ax)

#nmergedagain = mergedagain
v = ['id', 'Days after treatment', 'Type']
df = mergedagain.reset_index().set_index(v)._get_numeric_data().drop('index', axis=1)
df = df.loc[df.sum(axis=1) != 0, df.sum() != 0]
bresult = df.div(df.xs(0, level=1).droplevel(1), level=0).apply(np.log2)
bresult.replace([np.inf, -np.inf], np.nan, inplace=True)
bresult.fillna(0, inplace=True)
bmeanresult = bresult.droplevel(0).groupby(level=[0,1]).mean()
bbaseline = df.droplevel(2).xs(0, level=1)
bnotbaseline = df.droplevel(2).drop(0, level=1)

def bmandf(df, baseline):
    result = pd.Series()
    for i in df.columns:
        try:
            result[i] = mannwhitneyu(baseline[i],df[i])[1]
        except:
            result[i] = np.nan
    return result

pvals = bnotbaseline.groupby(level=1).apply(bmandf, (bbaseline))
#qvals = pd.DataFrame(fdrcorrection(pvals.values.flatten())[0].reshape(pvals.shape), index = pvals.index, columns = pvals.columns)
bsignificantMatrix = pvals.loc[:, (pvals < pval).any(axis=0)]
#significantMatrix = qvals.loc[:, (qvals < pval).any(axis=0)]
plotdf = bmeanresult[bsignificantMatrix.columns]
#sns.heatmap(plotdf.T, cmap='vlag', yticklabels=True)
#plotdf = bmeanresult.loc[bmeanresult.sum(axis=1) != 0, bmeanresult.sum() != 0]
plotdf = plotdf.loc[plotdf.sum(axis=1) != 0, plotdf.sum() != 0]
sns.heatmap(plotdf.T, cmap='vlag', yticklabels=True, center=0, square=True)

'''
