#!/usr/bin/env python
import functions
import pandas as pd
import matplotlib.pyplot as plt

meta = pd.read_csv("../data/meta.csv").set_index('Sample ID')
gutmsp = pd.read_csv("../data/gutmsp.csv", index_col=0).T
oralmsp = pd.read_csv("../data/oralmsp.csv", index_col=0).T
guttaxo = pd.read_csv("../data/guttaxo.csv", index_col=0)
oraltaxo = pd.read_csv("../data/oraltaxo.csv", index_col=0)

x = 'Type'
taxoType='genus'
variables = ['Days after treatment', 'Type', 'Patient_Id']

taxoOralmsp = oralmsp.T.join(oraltaxo[taxoType]).groupby(taxoType).sum().T
metaTaxoOralmsp = taxoOralmsp.join(meta[variables], how='inner').set_index(variables)
functions.relabund(metaTaxoOralmsp.sort_index())
plt.show()
'''
19, 20, 21, 23, 24, 21, 20, 19, 15, 9,7,6
'''

taxoGutmsp = gutmsp.T.join(guttaxo[taxoType]).groupby(taxoType).sum().T
metaTaxoGutmsp = taxoGutmsp.join(meta[variables], how='inner').set_index(variables)
functions.relabund(metaTaxoGutmsp.sort_index())
plt.show()

#meta[["id", "Type", "Days after treatment", "entero"]].to_csv("../results/metaentero")

cor, pval = functions.corr(metaTaxoGutmsp.add_prefix('gut '), metaTaxoOralmsp.add_prefix('oral '), min_unique=5)
thresh = 0.005
fcor, fpval = cor.loc[pval.lt(thresh).any(axis=1), pval.lt(thresh).any(axis=0)],pval.loc[pval.lt(thresh).any(axis=1), pval.lt(thresh).any(axis=0)]
functions.clustermap(fcor, fpval.lt(0.005))
plt.savefig('../results/oralgutcorrgenus.svg')
plt.show()

df1 = metaTaxoGutmsp.xs('DONOR', level=1).reset_index(drop=True) 
df2 = metaTaxoGutmsp.xs(0, level=0).droplevel(0)
df = pd.concat([df1,df2])
df = functions.mult(df)
df1["vals"] = 0
df2["vals"] = 1
grouping = pd.concat([df1,df2])['vals']

from skbio.stats.composition import ancom
from scipy.stats import mannwhitneyu
ancomdf, percentdf = ancom(
        df,
        grouping,
        significance_test=mannwhitneyu,
        multiple_comparisons_correction=None
)

import seaborn as sns
functions.heatmap(percentdf.loc[ancomdf["Reject null hypothesis"] == True])
plt.show()

'''
metaTaxoGutmsp.corr()[
    "Faecalibacterium prausnitzii 5"
].sort_values().dropna().to_frame().to_csv("fraus5.csv")
'''

edges = cor.stack().to_frame().reset_index().rename(columns={'level_0': 'source', 'level_1': 'target', 0: 'weight'})
edges = edges.loc[(edges.weight > 0.5) | (edges.weight < -0.5)]
edges.to_csv("tmp.csv")
G = nx.from_pandas_edgelist(edges)



# oral - no sif
df1 = metaTaxoOralmsp.xs(0, level=0).reset_index(drop=True) 
df2 = metaTaxoOralmsp.xs(7, level=0).droplevel(0)
df = pd.concat([df1,df2])
from skbio.stats.composition import multiplicative_replacement
df = pd.DataFrame(multiplicative_replacement(df), index=df.index, columns=df.columns)
df1["vals"] = 0
df2["vals"] = 1
grouping = pd.concat([df1,df2])['vals']
from skbio.stats.composition import ancom
from scipy.stats import mannwhitneyu
ancomdf, percentdf = ancom(df, grouping, significance_test=mannwhitneyu)
sns.heatmap(
    percentdf.loc[ancomdf["Reject null hypothesis"] == True], cmap="vlag", center=0
)
plt.show()

# gutbase - 
df1 = metaTaxoGutmsp.xs(0, level=0).reset_index(drop=True) 
df2 = metaTaxoGutmsp.xs([7,'FMT'], level=[0,1])
df = pd.concat([df1,df2])
from skbio.stats.composition import multiplicative_replacement
df = pd.DataFrame(multiplicative_replacement(df), index=df.index, columns=df.columns)
df1["vals"] = 0
df2["vals"] = 1
grouping = pd.concat([df1,df2])['vals']
from skbio.stats.composition import ancom
from scipy.stats import mannwhitneyu
ancomdf, percentdf = ancom(df, grouping, significance_test=mannwhitneyu)
sns.heatmap(
    percentdf.loc[ancomdf["Reject null hypothesis"] == True], cmap="vlag", center=0
)
plt.show()
