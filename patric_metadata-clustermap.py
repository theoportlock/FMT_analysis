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
from statsmodels.stats.multitest import multipletests

samples_metadata = pd.read_csv('metadata.csv', index_col=0).dropna()._get_numeric_data()
samples_patric = pd.read_csv("patricNorm.csv", index_col=0).T

# Join patric information
samples_patricMetadata = samples_patric.join(samples_metadata, how='inner')

# Calculate and format correlations
correlationArray, uncorrectedPValueArray = spearmanr(samples_patricMetadata)
correlations = pd.DataFrame(
    correlationArray,
    index=samples_patricMetadata.columns,
    columns=samples_patricMetadata.columns)
uncorrectedPValues = pd.DataFrame(
    uncorrectedPValueArray,
    index=samples_patricMetadata.columns,
    columns=samples_patricMetadata.columns)
slicedCorrelations = correlations.loc[
    samples_patric.columns,
    samples_metadata.columns]
slicedUncorrectedPValues = uncorrectedPValues.loc[
    samples_patric.columns,
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

#plt.show()
#plt.tight_layout()
plt.savefig("results/patric_metadata-clustermap.pdf")
