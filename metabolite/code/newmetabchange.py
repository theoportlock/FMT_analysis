#!/usr/bin/env python
from scipy import stats
from scipy.stats import mannwhitneyu
from scipy.stats import spearmanr
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

meta = pd.read_csv("../../data/newnewmetadata.csv").set_index('Sample ID')
df = pd.read_csv("../../data/metabolite.tsv", sep='\t')
#df = pd.read_csv("../../data/gutcard.csv", index_col=0).T
#df = pd.read_csv("../../data/gutcazy.csv", index_col=0).T
#df = pd.read_csv("../../data/gutkegg.csv", index_col=0).T
#df = pd.read_csv("../../data/gutpfam.csv", index_col=0).T
taxaType='species'
pval=0.05
#df = df.loc[:, df.nunique() > 20]
v=['Days after treatment', 'Type', 'id']
#df = df.join(meta[v], how='inner').set_index(v)

df['time_point'] = df['time_point'].str.replace('D', '')
df['time_point'] = df['time_point'].astype('float64') 
df.rename(columns={'time_point':'Days after treatment', 'IMP': 'Type'}, inplace=True)
df['Type'] = df['Type'].str.upper()
#df = df.set_index(v).join(meta.set_index(v)['MELD'], how='inner').set_index(v)
meta.id.replace('[a-zA-Z%]', '', regex=True, inplace=True)
meta.id = meta.id.astype('float64') 
df.id.replace('[a-zA-Z%]', '', regex=True, inplace=True)
df.id.replace('_', '', regex=True, inplace=True)
df.id = df.id.astype('float64') 
df = pd.merge(df, meta[v], on=v, how= 'inner').set_index(v)

# Delta B
bresult = df.div(df.xs(0, level=0).droplevel(0), level=2).apply(np.log2)
bresult.replace([np.inf, -np.inf], np.nan, inplace=True)
bresult.fillna(0, inplace=True)
bmeanresult = bresult.droplevel(2).groupby(level=[0,1]).mean().drop(0, level=0)
#sns.heatmap(bmeanresult.T,cmap='vlag', square=True)
bbaseline = df.droplevel(2).xs(0, level=0)
bnotbaseline = df.droplevel(2).drop(0, level=0)

def bmandf(df, baseline):
    result = pd.Series()
    for i in df.columns:
        try:
            result[i] = mannwhitneyu(baseline[i],df[i])[1]
        except:
            result[i] = np.nan
    return result

pvals = bnotbaseline.groupby(level=[0, 1]).apply(bmandf, (bbaseline))
#qvals = pd.DataFrame(fdrcorrection(pvals.values.flatten())[0].reshape(pvals.shape), index = pvals.index, columns = pvals.columns)
bsignificantMatrix = pvals.loc[:, (pvals < pval).any(axis=0)]
#significantMatrix = qvals.loc[:, (qvals < pval).any(axis=0)]
plotdf = bmeanresult[bsignificantMatrix.columns]
sns.heatmap(plotdf.sort_index(level=1).T, cmap='vlag', yticklabels=True, square=True)

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
#sns.heatmap(plotdf.sort_index(level=1).T, cmap='vlag', yticklabels=True, square=True)
sns.heatmap(pmeanresult.sort_index(level=1).T, cmap='vlag', yticklabels=True, square=True, center=0)

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
#qvals = pd.DataFrame(fdrcorrection(pvals.values.flatten())[0].reshape(pvals.shape), index = pvals.index, columns = pvals.columns)
dsignificantMatrix = pvals.loc[:, (pvals < pval).any(axis=0)]
#significantMatrix = qvals.loc[:, (qvals < pval).any(axis=0)]
plotdf = dmeanresult[dsignificantMatrix.columns]
sns.heatmap(plotdf.sort_index(level=1).T, cmap='vlag', yticklabels=True, square=True, center=0)

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
