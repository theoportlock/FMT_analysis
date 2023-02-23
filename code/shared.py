#!/usr/bin/env python
from scipy import stats
from scipy.stats import wilcoxon
from scipy.stats import spearmanr
from statsmodels.stats.multitest import fdrcorrection
from statannot import add_stat_annotation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

pval=0.1
taxaType='species'
meta = pd.read_csv("../../data/newnewmetadata.csv").set_index('Sample ID')
gmsp = pd.read_csv("../../data/gutmsp.csv", index_col=0)
omsp = pd.read_csv("../../data/oralmsp.csv", index_col=0)
gtaxo = pd.read_csv("../../data/guttaxo.csv", index_col=0)
otaxo = pd.read_csv("../../data/oraltaxo.csv", index_col=0)

gmsptaxo = gmsp.join(gtaxo[taxaType], how='inner').groupby(taxaType).sum().T
omsptaxo = omsp.join(otaxo[taxaType], how='inner').groupby(taxaType).sum().T
gmsptaxo = gmsptaxo[gmsptaxo.columns[~gmsptaxo.columns.str.contains('unclassified')]]
omsptaxo = omsptaxo[omsptaxo.columns[~omsptaxo.columns.str.contains('unclassified')]]
gmsptaxo.columns = gmsptaxo.columns.str.replace(r' /.*', '', regex=True)
omsptaxo.columns = omsptaxo.columns.str.replace(r' /.*', '', regex=True)
gmsptaxo.columns = gmsptaxo.columns.str.replace(r' &.*', '', regex=True)
omsptaxo.columns = omsptaxo.columns.str.replace(r' &.*', '', regex=True)
gmsptaxo = gmsptaxo.loc[:, gmsptaxo.nunique() > 20]
omsptaxo = omsptaxo.loc[:, omsptaxo.nunique() > 20]
v=['Days after treatment', 'Type', 'id']
gdf = gmsptaxo.join(meta[v], how='inner').set_index(v)
odf = omsptaxo.join(meta[v], how='inner').set_index(v)
gvals = gdf.reset_index().groupby('id').count()['Days after treatment']
ovals = odf.reset_index().groupby('id').count()['Days after treatment']
gfvals = gvals.loc[(gvals == 4) | (gvals == 0)].index
ofvals = ovals.loc[(ovals == 4) | (ovals == 0)].index
vals = set.intersection(set(gfvals), set(ofvals))
gdf.reset_index(inplace=True)
odf.reset_index(inplace=True)
gdf = gdf[gdf['id'].isin(vals)].set_index(v)
odf = odf[odf['id'].isin(vals)].set_index(v)
gdf = gdf.loc[:, gdf.nunique() > 20].T
odf = odf.loc[:, odf.nunique() > 20].T
sgdf = gdf.copy()
sodf = odf.copy()
sgdf['Site'] = 'gut'
sodf['Site'] = 'oral'
sgdf.set_index('Site', append=True, inplace=True)
sodf.set_index('Site', append=True, inplace=True)
gdf = gdf.T.add_prefix('gut ')
odf = odf.T.add_prefix('oral ')

df = sgdf.T.join(sodf.T)
#df.columns = df.columns.to_flat_index()
jdf = gdf.join(odf, how='inner')


# Delta B
bresult = jdf.div(jdf.xs(0, level=0).droplevel(0), level=2).apply(np.log2)
bresult.replace([np.inf, -np.inf], np.nan, inplace=True)
bresult.fillna(0, inplace=True)
bmeanresult = bresult.droplevel(2).groupby(level=[0,1]).mean().drop(0, level=0)
bbaseline = jdf.droplevel(2).xs(0, level=0)
bnotbaseline = jdf.droplevel(2).drop(0, level=0)

def bmandf(df, baseline):
    result = pd.Series()
    baselinevar = df.reset_index()['Type'].unique()[0]
    for i in baseline.columns:
        try:
            result.loc[i] = mannwhitneyu(baseline.loc[baselinevar,i],df[i])[1]
        except:
            result.loc[i] = np.nan
    return result

pvals = bnotbaseline.groupby(level=[0, 1]).apply(bmandf, (bbaseline))
#qvals = pd.DataFrame(fdrcorrection(pvals.values.flatten())[0].reshape(pvals.shape), index = pvals.index, columns = pvals.columns)
bsignificantMatrix = pvals.loc[:, (pvals < pval).any(axis=0)]
#significantMatrix = qvals.loc[:, (qvals < pval).any(axis=0)]
plotdf = bmeanresult[bsignificantMatrix.columns]
sns.heatmap(plotdf.sort_index(level=1).T.sort_index(), cmap='vlag', yticklabels=True, square=True, center=0)

# Delta P
presult = jdf.xs('FMT',level=1).droplevel(1).groupby(level=0).mean().div(
    jdf.xs('PLACEBO',level=1).droplevel(1).groupby(level=0).mean()).apply(np.log2)
pmeanresult = presult.replace([np.inf, -np.inf], np.nan).fillna(0)
pbaseline = jdf.xs('FMT', level=1).droplevel(1)
pnotbaseline = jdf.xs('PLACEBO', level=1).droplevel(1)

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
sns.heatmap(plotdf.sort_index(level=1).T.sort_index(), cmap='vlag', yticklabels=True, square=True, center=0)

# Delta S
sresult = jdf.xs('FMT',level=1).droplevel(1).groupby(level=0).mean().div(
    jdf.xs('PLACEBO',level=1).droplevel(1).groupby(level=0).mean()).apply(np.log2)
smeanresult = sresult.replace([np.inf, -np.inf], np.nan).fillna(0)
sbaseline = jdf.xs('FMT', level=1).droplevel(1)
snotbaseline = jdf.xs('PLACEBO', level=1).droplevel(1)

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

pvals = pd.DataFrame(snotbaseline.groupby(level=0).apply(smandf, (pbaseline)))
#qvals = pd.DataFrame(fdrcorrection(pvals.values.flatten())[0].reshape(pvals.shape), index = pvals.index, columns = pvals.columns)
ssignificantMatrix = pvals.loc[:, (pvals < pval).any(axis=0)]
#significantMatrix = qvals.loc[:, (qvals < pval).any(axis=0)]
plotdf = smeanresult[ssignificantMatrix.columns]
sns.heatmap(plotdf.sort_index(level=1).T.sort_index(), cmap='vlag', yticklabels=True, square=True, center=0)

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


'''
#corheatmap
cordf = df.corr(method='spearman')
fcordf = cordf.loc[np.abs(cordf > 0.75).any(),np.abs(cordf > 0.75).any()]
sns.clustermap(fcordf, xticklabels=True, yticklabels=True, cmap='vlag', center=0)
plt.show()
'''
'''
#shared
staxo = gtaxo.loc[gtaxo.species.isin(otaxo.species)]
gmsp.join
gmsptaxo = gmsp.join(gtaxo[taxaType], how='inner').groupby(taxaType).sum()
shared = gmsptaxo.loc[gmsptaxo.index.isin(otaxo.species)].sum()
allgut = gmsptaxo.sum()
shared = shared / allgut
meta['shared'] = shared
meta.rename(columns={'Type': 'Arm'}, inplace=True)
df = meta.loc[meta.Arm != 'DONOR']
df = df.dropna(subset=['shared'])

x='Days after treatment'
y='shared'
hue='Arm'
plt.rcParams["figure.figsize"] = (4,4)
fig, ax = plt.subplots(1,1)
sns.boxplot(data=df, x=x, y=y, hue=hue, showfliers=False, linewidth=1, ax=ax)
sns.stripplot(data=df, x='Days after treatment', y='shared', hue='Arm', dodge=True,size=3, linewidth=1, ax=ax)

# Stats
import itertools
box_pairs = list(itertools.combinations([(a,b) for a in df[x].unique() for b in df[hue].unique()],2))
test_short_name = 'M.W.W'
pvalues = {}
for pair in box_pairs:
    data1 = df.groupby([x,hue])[y].get_group(pair[0])
    data2 = df.groupby([x,hue])[y].get_group(pair[1])
    stat, p = mannwhitneyu(data1, data2)
    pvalues[pair] = p
pvals = pd.Series(pvalues)
#qvals = pd.Series(fdrcorrection(pvals)[1], index=pvals.index)
#sigqvals = qvals.loc[pvals < 0.06]
sigqvals = pvals.loc[pvals < 0.06]

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
#df.groupby(['Days after treatment', 'Arm']).mean().sort_index(level=1).shared.plot.bar()

cordf = df.corr(method='spearman')
for i in cordf.columns:
    cordf.loc[i,i] = 0
fcordf = cordf.loc[np.abs(cordf > 0.5).any(),np.abs(cordf > 0.5).any()]
sns.clustermap(fcordf, xticklabels=True, yticklabels=True, cmap='vlag', center=0)

plt.show()
df.replace([np.inf, -np.inf], np.nan, inplace=True)
df.dropna(subset=['MELD'], inplace=True)
df.dropna(axis=1,inplace=True)
