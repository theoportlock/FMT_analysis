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
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

samples_metadata = pd.read_csv('metadata.csv', index_col=0).dropna()
samples_patric = pd.read_csv("patricNorm.csv", index_col=0).T

# Join patric information
samples_patricMetadata = samples_patric.join(samples_metadata['Days after treatment'], how='inner')

# Calculate and format changes
averageDays = samples_patricMetadata.groupby('Days after treatment').mean()
averageDays = averageDays.loc[(averageDays.sum(axis=1) != 0), (averageDays.sum(axis=0) != 0)]

baseline = samples_patricMetadata.loc[samples_patricMetadata['Days after treatment'] == 0]
difference = averageDays.div(baseline.mean())
difference = difference.loc[(difference.sum(axis=1) != 0), (difference.sum(axis=0) != 0)]
difference = difference.loc[(difference.sum(axis=1) != 1), (difference.sum(axis=0) != 1)]

stats = samples_patricMetadata.groupby('Days after treatment').apply(lambda group: mannwhitneyu(baseline, group)[1])
statsdf = pd.DataFrame(stats.to_list(), columns=samples_patricMetadata.columns, index=averageDays.index).drop(['Days after treatment'], axis=1)
statsdf.loc[:, (statsdf < 0.05).any(axis=0)]




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
