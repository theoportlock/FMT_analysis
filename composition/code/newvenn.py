#!/usr/bin/env python
from scipy import stats
from scipy.stats import wilcoxon
from scipy.stats import spearmanr
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

meta = pd.read_csv("../../data/newnewmetadata.csv").set_index('Sample ID')
msp = pd.read_csv("../../data/gutmsp.csv", index_col=0)
taxaType='genus'
pval=0.05
taxo = pd.read_csv("../../data/guttaxo.csv", index_col=0)
msptaxo = msp.join(taxo[taxaType], how='inner').groupby(taxaType).sum().T
msptaxo = msptaxo[msptaxo.columns[~msptaxo.columns.str.contains('unclassified')]]
msptaxo.columns = msptaxo.columns.str.replace(r' /.*', '', regex=True)
msptaxo.columns = msptaxo.columns.str.replace(r' &.*', '', regex=True)
msptaxo = msptaxo.loc[:, msptaxo.nunique() > 20]
v=['Days after treatment', 'Type', 'id']
df = msptaxo.join(meta[v], how='inner').set_index(v)
#only look for 4 timepoints
vals = df.reset_index().groupby('id').count()['Days after treatment']
filteredvals = vals.loc[(vals == 4) | (vals == 0)].index
df.reset_index(inplace=True)
df = df[df['id'].isin(filteredvals)].set_index(v)
df = df.loc[:, df.nunique() > 20]

# Delta B
bresult = df.drop('DONOR', level=1).div(df.xs(0, level=0).droplevel(0), level=2).apply(np.log2)
bresult.replace([np.inf, -np.inf], np.nan, inplace=True)
bresult.fillna(0, inplace=True)
bmeanresult = bresult.droplevel(2).groupby(level=[0,1]).mean().drop(0, level=0)
bbaseline = df.droplevel(2).xs(0, level=0)
bnotbaseline = df.droplevel(2).drop(0, level=0).drop('DONOR', level=1)

def bmandf(df, baseline):
    result = pd.Series()
    baselinevar = df.reset_index()['Type'].unique()[0]
    for i in df.columns:
        try:
            result[i] = wilcoxon(baseline.loc[baselinevar,i],df[i])[1]
        except:
            result[i] = np.nan
    return result

pvals = bnotbaseline.groupby(level=[0, 1]).apply(bmandf, (bbaseline))
#qvals = pd.DataFrame(fdrcorrection(pvals.values.flatten())[0].reshape(pvals.shape), index = pvals.index, columns = pvals.columns)
bsignificantMatrix = pvals.loc[:, (pvals < pval).any(axis=0)]
#significantMatrix = qvals.loc[:, (qvals < pval).any(axis=0)]
plotdf = bmeanresult[bsignificantMatrix.columns]
sns.heatmap(plotdf.sort_index(level=1).T.sort_index(), cmap='vlag', yticklabels=True, square=True)

# Delta P
presult = df.xs('FMT',level=1).droplevel(1).groupby(level=0).mean().div(
    df.xs('PLACEBO',level=1).droplevel(1).groupby(level=0).mean()).apply(np.log2)
pmeanresult = presult.replace([np.inf, -np.inf], np.nan).fillna(0)
pbaseline = df.xs('FMT', level=1).droplevel(1)
pnotbaseline = df.xs('PLACEBO', level=1).droplevel(1)

def pmandf(df, pbaseline):
    var = df.index.unique()[0]
    base = pbaseline.loc[var]
    result = pd.Series()
    for i in base.columns:
        try:
            result[i] = mannwhitneyu(base[i],df[i])[1]
        except:
            result[i] = np.nan
    return result

pvals = pd.DataFrame(pnotbaseline.groupby(level=0).apply(pmandf, (pbaseline)))
#qvals = pd.DataFrame(fdrcorrection(pvals.values.flatten())[0].reshape(pvals.shape), index = pvals.index, columns = pvals.columns)
psignificantMatrix = pvals.loc[:, (pvals < pval).any(axis=0)]
#significantMatrix = qvals.loc[:, (qvals < pval).any(axis=0)]
plotdf = pmeanresult[psignificantMatrix.columns]
sns.heatmap(plotdf.sort_index(level=1).T.sort_index(), cmap='vlag', yticklabels=True, square=True)

# Delta D
dresult = df.droplevel(2).groupby(level=[0,1]).mean().div(
    df.droplevel(2).xs('DONOR', level=1).mean()).apply(np.log2)
dmeanresult = dresult.replace([np.inf, -np.inf], np.nan).fillna(0)
dbaseline = df.droplevel(2).xs('DONOR', level=1)
dnotbaseline = df.droplevel(2).drop('DONOR', level=1)
def dmandf(df, dbaseline):
    var = df.index.unique()[0]
    result = pd.Series()
    for i in df.columns:
        try:
            result[i] = mannwhitneyu(dbaseline[i],df[i])[1]
        except:
            result[i] = np.nan
    return result

pvals = pd.DataFrame(dnotbaseline.groupby(level=[0,1]).apply(dmandf, (dbaseline)))
# for just 0 timepoint
pvals = pvals.xs(0, level=0)
#qvals = pd.DataFrame(fdrcorrection(pvals.values.flatten())[0].reshape(pvals.shape), index = pvals.index, columns = pvals.columns)
dsignificantMatrix = pvals.loc[:, (pvals < pval).any(axis=0)]
#significantMatrix = qvals.loc[:, (qvals < pval).any(axis=0)]
plotdf = dmeanresult[dsignificantMatrix.columns]
sns.heatmap(plotdf.sort_index(level=1).T.sort_index(), cmap='vlag', yticklabels=True, square=True)
'''
sns.stripplot(data=pmeanresult[inter].T.sort_index().T.reset_index().melt(id_vars='Days 
after treatment'), x='variable', y='value', hue='Days after treatment')
plt.xticks(rotation=90)
plt.show()
'''

# intersection
inter = set(dsignificantMatrix.columns) & set(bsignificantMatrix) & set(psignificantMatrix)
sns.heatmap(dmeanresult[inter].sort_index(level=1).T.sort_index(), cmap='vlag', yticklabels=True, xticklabels=True, center=0, square=True)
#plt.savefig('d.svg')
plt.show()
sns.heatmap(pmeanresult[inter].T.sort_index(), cmap='vlag', yticklabels=True, center=0, square=True)
#plt.savefig('p.svg')
plt.show()
sns.heatmap(bmeanresult[inter].sort_index(level=1).T.sort_index(), cmap='vlag', yticklabels=True, center=0, square=True)
#plt.savefig('b.svg')
plt.show()

# Venn
from matplotlib_venn import venn3
from itertools import combinations
def combs(x): return [c for i in range(1, len(x)+1) for c in combinations(x,i)]
comb = combs([dsignificantMatrix, bsignificantMatrix, psignificantMatrix])
result = []
for i in comb:
    if len(i) > 1:
        result.append(len(set.intersection(*(set(j.columns) for j in i))))
    else:
        result.append(len(i[0].columns))

venn3(subsets = result)
plt.show()
