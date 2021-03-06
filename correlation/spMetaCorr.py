#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Plotting an annotated seaborn correlation heatmap between MSP data and metadata
Theo Portlock
%autoindent
'''

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

sns.set(rc={'figure.figsize':(11.7,8.27)})
samples_metadata = pd.read_csv('../data/cleanMeta.csv').set_index('Sample ID').fillna(0)._get_numeric_data()
othermeta = pd.read_csv('../data/plasmaBileAcid.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/plasmaTryptophan.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/stoolBileAcid.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/stoolNmr.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/stoolTryptophan.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/urineBileAcid.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/urineNmr.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/urineTryptophan.tsv', sep= '\t').set_index('Sample ID')

#samples_metadata.drop("Creatinine", axis=1, inplace=True)
samples_metadata = samples_metadata.join(othermeta)
samples_metadata = samples_metadata.loc[samples_metadata.MELDNA != 0]
samples_metadata.dropna(inplace=True)
#msp_taxonomy = pd.read_csv("../data/guttaxo.csv", index_col=0)
#msp_samples = pd.read_csv("../data/gutmsp.csv", index_col=0)
msp_taxonomy = pd.read_csv("../data/oraltaxo.csv", index_col=0)
msp_samples = pd.read_csv("../data/oralmsp.csv", index_col=0)

# Join taxo information
taxaType='species'
taxonomy_samples = msp_samples.join(msp_taxonomy[taxaType], how='inner').groupby(taxaType).sum()
samples_taxonomy = taxonomy_samples.T
samples_taxonomyMetadata = samples_taxonomy.join(samples_metadata, how='inner').astype('float')

# Calculate and format correlations
correlationArray, uncorrectedPValueArray = spearmanr(samples_taxonomyMetadata)
correlations = pd.DataFrame(
    correlationArray,
    index=samples_taxonomyMetadata.columns,
    columns=samples_taxonomyMetadata.columns)
uncorrectedPValues = pd.DataFrame(
    uncorrectedPValueArray,
    index=samples_taxonomyMetadata.columns,
    columns=samples_taxonomyMetadata.columns)
slicedCorrelations = correlations.loc[
    samples_taxonomy.columns,
    samples_metadata.columns]
slicedUncorrectedPValues = uncorrectedPValues.loc[
    samples_taxonomy.columns,
    samples_metadata.columns]
filteredSlicedUncorrectedPValues = slicedUncorrectedPValues.loc[
    (slicedCorrelations.sum(axis=1) != 0),
    (slicedCorrelations.sum(axis=0) != 0)]
filteredSlicedCorrelations = slicedCorrelations.loc[
    (slicedCorrelations.sum(axis=1) != 0),
    (slicedCorrelations.sum(axis=0) != 0)]
significantMatrix = pd.DataFrame(
    multipletests(filteredSlicedUncorrectedPValues.values.flatten())[0].reshape(filteredSlicedUncorrectedPValues.shape),
    index = filteredSlicedUncorrectedPValues.index,
    columns = filteredSlicedUncorrectedPValues.columns)

# Plot
g = sns.clustermap(
    #filteredSlicedCorrelations,
    filteredSlicedCorrelations.loc[significantMatrix.any(axis=1),:],
    cmap="vlag",
    vmin=-1,
    vmax=1, 
    yticklabels=True,
    xticklabels=True)

for tick in g.ax_heatmap.get_yticklabels():
    tick.set_rotation(0)

for i, ix in enumerate(g.dendrogram_row.reordered_ind):
    for j, jx in enumerate(g.dendrogram_col.reordered_ind):
        text = g.ax_heatmap.text(
            j + 0.5,
            i + 0.5,
            #"*" if significantMatrix.values[ix, jx] else "",
            "*" if significantMatrix.loc[significantMatrix.any(axis=1),:].values[ix, jx] else "",
            ha="center",
            va="center",
            color="black",
        )
        text.set_fontsize(8)

g.ax_heatmap.tick_params(axis='both', labelsize=8)
plt.show()
#plt.tight_layout()
#plt.savefig("../results/" + taxaType + "_clustermap.pdf")
