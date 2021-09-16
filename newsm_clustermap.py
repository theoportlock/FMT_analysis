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

#samples_metadata = pd.read_csv('metadata.csv', index_col=0).dropna()._get_numeric_data()
#samples_metadata = pd.read_csv('newmergedmetadata.csv').drop_duplicates(subset='ID_x', keep='first').set_index('ID_x').sort_index().dropna(axis=1)._get_numeric_data()
#samples_metadata = pd.read_csv('newmergedmetadata.csv').drop_duplicates(subset='ID_x', keep='first').set_index('ID_x').sort_index().fillna(0)._get_numeric_data()
#samples_metadata = pd.read_csv('newmergedmetadata.csv').drop_duplicates(subset='ID_x', keep='first').set_index('ID_x').sort_index()._get_numeric_data()
samples_metadata = pd.read_csv('newnewmetadata.csv').drop_duplicates(subset='Sample ID', keep='first').set_index('Sample ID').sort_index().dropna(thresh=20, axis=1)._get_numeric_data()
samples_sm = pd.read_csv("antismashNorm.csv", index_col=0).T
samples_sm.columns.name = 'Secondary metabolite'

# Join sm information
samples_smMetadata = samples_sm.join(samples_metadata, how='inner')

# Calculate and format correlations
correlationArray, uncorrectedPValueArray = spearmanr(samples_smMetadata, nan_policy='omit')
correlations = pd.DataFrame(
    correlationArray,
    index=samples_smMetadata.columns,
    columns=samples_smMetadata.columns)
uncorrectedPValues = pd.DataFrame(
    uncorrectedPValueArray,
    index=samples_smMetadata.columns,
    columns=samples_smMetadata.columns)
slicedCorrelations = correlations.loc[
    samples_sm.columns,
    samples_metadata.columns]
slicedUncorrectedPValues = uncorrectedPValues.loc[
    samples_sm.columns,
    samples_metadata.columns]
filteredSlicedUncorrectedPValues = slicedUncorrectedPValues.loc[
    (slicedCorrelations.sum(axis=1) != 0),
    (slicedCorrelations.sum(axis=0) != 0)]
filteredSlicedCorrelations = slicedCorrelations.loc[
    (slicedCorrelations.sum(axis=1) != 0),
    (slicedCorrelations.sum(axis=0) != 0)]
significantMatrix = pd.DataFrame(
    fdrcorrection(filteredSlicedUncorrectedPValues.values.flatten())[0].reshape(filteredSlicedUncorrectedPValues.shape),
    index = filteredSlicedUncorrectedPValues.index,
    columns = filteredSlicedUncorrectedPValues.columns)

# Plot
g = sns.clustermap(
    filteredSlicedCorrelations,
    #filteredSlicedCorrelations.loc[significantMatrix.any(axis=1),:],
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
            "*" if significantMatrix.values[ix, jx] else "",
            #"*" if significantMatrix.loc[significantMatrix.any(axis=1),:].values[ix, jx] else "",
            ha="center",
            va="center",
            color="black",
        )
        text.set_fontsize(8)

#plt.show()
#plt.tight_layout()
plt.savefig("results/sm_metadata-clustermap.pdf")
