#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Plotting an annotated seaborn correlation heatmap between MSP data and metadata
Theo Portlock
'''
%autoindent

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import fdrcorrection

samples_metadata = pd.read_csv("newmergedmetadata.csv").set_index('ID_x')
gut_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
gut_msp_data = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
oral_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)
oral_msp_data = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)

# Format data
taxa_type='genus'

# donor
dtaxa_not_sum = gut_msp_data.join(gut_taxonomy[taxa_type], how='inner').set_index(taxa_type).T
otaxa_not_sum = oral_msp_data.join(oral_taxonomy[taxa_type], how='inner').set_index(taxa_type).T
otaxa = otaxa_not_sum.groupby(otaxa_not_sum.columns, axis=1).sum()
dtaxa = dtaxa_not_sum.groupby(dtaxa_not_sum.columns, axis=1).sum()
sharedind = pd.Index.intersection(otaxa.index, dtaxa.index)
otaxa, dtaxa = otaxa.loc[sharedind], dtaxa.loc[sharedind]
otaxa.columns.name, dtaxa.columns.name = 'Oral', 'Gut'

#correlationArray, uncorrectedPValues = spearmanr(otaxa.loc[sharedind], dtaxa.loc[sharedind], axis=0)
correlationArray, uncorrectedPValueArray = spearmanr(otaxa, dtaxa, axis=0)
correlations = pd.DataFrame(
    correlationArray,
    index=otaxa.columns.append(dtaxa.columns),
    columns=otaxa.columns.append(dtaxa.columns))
uncorrectedPValues = pd.DataFrame(
    uncorrectedPValueArray,
    index=otaxa.columns.append(dtaxa.columns),
    columns=otaxa.columns.append(dtaxa.columns))

correlations = correlations.dropna(thresh=3,axis=1).dropna(thresh=3, axis=0)
uncorrectedPValues = uncorrectedPValues.dropna(thresh=3,axis=1).dropna(thresh=3, axis=0)

slicedCorrelations = correlations.iloc[
        len(otaxa.columns):,
        :len(otaxa.columns)]
slicedUncorrectedPValues = uncorrectedPValues.iloc[
        len(otaxa.columns):,
        :len(otaxa.columns)]

significantMatrix = pd.DataFrame(
    fdrcorrection(slicedUncorrectedPValues.values.flatten())[0].reshape(slicedUncorrectedPValues.shape),
    index = slicedUncorrectedPValues.index,
    columns = slicedUncorrectedPValues.columns)

# Plot
#sns.set()
#sns.set_theme()
sns.set(rc={'figure.figsize':(10,10)})
#sns.set_context("paper", rc={"font.size":3,"axes.titlesize":3,"axes.labelsize":10})
sns.set_theme(font_scale=0.3)

g = sns.clustermap(
    slicedCorrelations,
    #slicedCorrelations.loc[significantMatrix.any(axis=1),:],
    cmap="vlag",
    vmin=-1,
    vmax=1, 
    #cbar_kws=dict(use_gridspec=False,shrink=0.25),
    #cbar=False,
    yticklabels=True,
    xticklabels=True)

for tick in g.ax_heatmap.get_yticklabels():
    tick.set_rotation(0)

for i, ix in enumerate(g.dendrogram_row.reordered_ind):
    for j, jx in enumerate(g.dendrogram_col.reordered_ind):
        text = g.ax_heatmap.text(
            j + 0.5,
            i + 0.5,
            "*" if significantMatrix.values[ix, jx] else "",
            #"*" if significantMatrix.loc[significantMatrix.any(axis=1),:].values[ix, jx] else "",
            ha="center",
            va="center",
            color="black",
        )
        text.set_fontsize(3)

g.ax_heatmap.set(xlabel="Oral genus", ylabel="Gut genus")
plt.subplots_adjust(top=0.98, left=0.02, right=0.9, bottom=0.16)
g.ax_cbar.set_position([0.001, 0.85, g.ax_row_dendrogram.get_position().width, 0.05])
g.ax_cbar.set_title('FDR adjusted\nSpearman Correlation')
g.ax_cbar.tick_params(axis='y', length=5)
#plt.subplots_adjust(right=0.9, bottom=0.16)
#g.ax_cbar.autoscale_view()

#plt.show()
#plt.tight_layout(h_pad=0.1, w_pad=0.1)
#plt.savefig("results/oral_gut_clustermap.pdf")


'''
#network stuff
edges = slicedCorrelations.stack().reset_index().set_index(['level_0','level_1'])
sigvals = significantMatrix.stack().reset_index().set_index(['level_0','level_1'])
sigedges = edges[sigvals].dropna()
#gx = nx.from_pandas_edgelist(sigedges.reset_index(), 'level_0', 'level_1', 0)
#sigedges.to_csv('edges.csv')
sigedges = sigedges.reset_index()
sigedges['source'] = sigedges['level_0'].apply(lambda x: f"oral {x}")
sigedges['target'] = sigedges['level_1'].apply(lambda x: f"gut {x}")
sigedges[['source', 'target', 0]].to_csv('edges.csv')
df1 = sigedges[['source','level_0']].drop_duplicates().set_axis(['name1', 'name2'], axis=1)
df2 = sigedges[['target', 'level_1']].drop_duplicates().set_axis(['name1', 'name2'], axis=1)
df1['Site'], df2['Site'] = 'Oral', 'Gut'
renamedf = pd.concat([df1, df2])
renamedf.to_csv('renamedf.csv')
'''
