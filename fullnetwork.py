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

taxaType='species'
samples_metadata = pd.read_csv("newnewmetadata.csv").set_index('Sample ID').dropna(axis=0, subset=['MELD'])

gmsp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
gmsp_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
samples_gtaxonomy = gmsp_samples.join(gmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
gmspdf = samples_gtaxonomy.add_prefix('gut ')

omsp_samples = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)
omsp_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)
samples_otaxonomy = omsp_samples.join(omsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
omspdf = samples_otaxonomy.add_prefix('oral ')

df = omspdf.join(gmspdf, how='inner')

correlationArray, uncorrectedPValueArray = spearmanr(df)
correlations = pd.DataFrame(
    correlationArray,
    index=df.columns,
    columns=df.columns)
uncorrectedPValues = pd.DataFrame(
    uncorrectedPValueArray,
    index=df.columns,
    columns=df.columns)

correlations = correlations.dropna(thresh=3,axis=1).dropna(thresh=3, axis=0)
uncorrectedPValues = uncorrectedPValues.dropna(thresh=3,axis=1).dropna(thresh=3, axis=0)

slicedCorrelations = correlations
slicedUncorrectedPValues = uncorrectedPValues

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
slicedCorrelations.columns.name='target1'
significantMatrix.columns.name='target1'
edges = slicedCorrelations.stack().reset_index().set_index(['species','target1'])
sigvals = significantMatrix.stack().reset_index().set_index(['species','target1'])
sigedges = edges[sigvals].dropna()
import networkx as nx
gx = nx.from_pandas_edgelist(sigedges.reset_index(), 'species', 'target1', 0)
#sigedges.to_csv('fulledges.csv')
sigedges = sigedges.reset_index()
sigedges['source'] = sigedges['species'].apply(lambda x: f"gut {x}")
sigedges['target'] = sigedges['target1'].apply(lambda x: f"oral {x}")
#sigedges[['source', 'target', 0]].to_csv('fulledges.csv')
df1 = sigedges[['source','species']].drop_duplicates().set_axis(['name1', 'name2'], axis=1)
df2 = sigedges[['target', 'target1']].drop_duplicates().set_axis(['name1', 'name2'], axis=1)
df1['Site'], df2['Site'] = 'Oral', 'Gut'
renamedf = pd.concat([df1, df2])
renamedf.to_csv('fullrenamedf.csv')
'''
