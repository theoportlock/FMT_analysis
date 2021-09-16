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
sns.set_theme(font_scale=0.8)
'''
Cant use significance for the placebo group as the group is too small
'''

samples_metadata = pd.read_csv('metadata.csv', index_col=0)
card_samples = pd.read_csv("cazymeNorm.csv", index_col=0)
#var = 'AMR Gene Family'
samples_card = card_samples[[*samples_metadata.loc[samples_metadata['Type'] == 'FMT'].index]].T
plac_samples_card = card_samples[[*samples_metadata.loc[samples_metadata['Type'] == 'PLACEBO'].index]].T
#samples_card = card_samples[[*samples_metadata.loc[samples_metadata['Type'] == 'PLACEBO'].index, var]].groupby(var).sum().T

# Join card information
samples_cardMetadata = samples_card.join(samples_metadata['Days after treatment'], how='inner')

# Calculate and format changes
averageDays = samples_cardMetadata.groupby('Days after treatment').mean()
#averageDays = np.arctan(averageDays.loc[(averageDays.sum(axis=1) != 0), (averageDays.sum(axis=0) != 0)])
#averageDays = averageDays.loc[(averageDays.sum(axis=1) != 0), (averageDays.sum(axis=0) != 0)]
#averageDays = np.log1p(averageDays)

baseline = samples_cardMetadata.loc[samples_cardMetadata['Days after treatment'] == 0]
#difference = averageDays.div(baseline.mean())
difference = averageDays.sub(baseline.mean())
#difference = difference.loc[(difference.sum(axis=1) != 0), (difference.sum(axis=0) != 0)]
#difference = difference.loc[(difference.sum(axis=1) != 1), (difference.sum(axis=0) != 1)]
#difference.replace([np.inf, -np.inf], np.nan, inplace=True)
#difference.dropna(axis=1, inplace=True)
#difference = np.tanh(difference * 1e4)

def mandf(df):
    base = baseline.loc[(baseline.sum(axis=1)!=0),(baseline.sum(axis=0)!=0)].copy()
    day = df.loc[(df.sum(axis=1)!=0),(df.sum(axis=0)!=0)].copy()
    result = pd.Series()
    for i in base.columns:
        try:
            result[i] = mannwhitneyu(base[i],day[i])[1]
        except:
            result[i] = 1
    return result

stats = samples_cardMetadata.groupby('Days after treatment').apply(mandf)
significantMatrix = stats.loc[:, (stats < 0.01).any(axis=0)]
#plotMatrix = averageDays[significantMatrix.columns]
plotMatrix = difference[significantMatrix.columns]
plotMatrix = plotMatrix.T.iloc[:,1:].sort_values(7)
significantMatrix = significantMatrix.T.iloc[:,1:]
significantMatrix = significantMatrix.loc[plotMatrix.index]

# Plot
g = sns.clustermap(
#g = sns.heatmap(
    #data=difference.loc[:,significantMatrix].T.iloc[:,1:],
    #data=difference.T,
    #data=averageDays.T,
    #data=difference.loc[:, stats.loc[, (stats <= 0.05).any(axis=0)].columns].T,
    #data=difference.loc[:, stats.loc[:, (stats < 0.05).any(axis=0)].columns].T,
    #data=averageDays.loc[:, stats.loc[:, (stats <= 0.05).any(axis=0)].columns].T,
    #data=plotMatrix.T.iloc[:,1:].sort_values(7),
    data=plotMatrix,
    #data=averageDays.loc[:,significantMatrix.columns].T,
    col_cluster=False,
    row_cluster=False,
    cmap="vlag",
    center=0,
    #z_score=0,
    #metric='correlation',
    #standard_scale=0,
    #method='single',
    #robust=True,
    #cmap="Blues",
    #vmin=-1,
    #vmax=1, 
    yticklabels=True,
    xticklabels=True)

for tick in g.ax_heatmap.get_yticklabels():
    tick.set_rotation(0)

annot=pd.DataFrame()
annot[(significantMatrix < 0.0005) & (plotMatrix > 0)] = '+++'
annot[(significantMatrix < 0.005) & (plotMatrix > 0)] = '++'
annot[(significantMatrix < 0.05) & (plotMatrix > 0)] = '+'
annot[(significantMatrix < 0.0005) & (plotMatrix < 0)] = '---'
annot[(significantMatrix < 0.005) & (plotMatrix < 0)] = '--'
annot[(significantMatrix < 0.05) & (plotMatrix < 0)] = '-'
annot[significantMatrix >= 0.05] = ''

#for i, ix in enumerate(g.dendrogram_row.reordered_ind):
#for j, ix in enumerate(list(range(len(plotMatrix.columns)))):
#    for j, jx in enumerate(list(range(len(plotMatrix.index)))):
for i, ix in enumerate(plotMatrix.index):
    for j, jx in enumerate(plotMatrix.columns):
        text = g.ax_heatmap.text(
        #text = g.text(
            j + 0.5,
            i + 0.5,
            #"*" if significantMatrix.values[i,j] <= 0.05 else "",
            #"+" if significantMatrix.values[jx,ix] <= 0.05 and plotMatrix.values[jx,ix] > 0 else "",
            annot.values[i,j],
            ha="center",
            va="center",
            color="black",
        )
        text.set_fontsize(8)

plt.show()
#plt.tight_layout()
#plt.savefig("results/cazymeHeatSigFmt.pdf")
 lt.savefig("results/cazymeHeatSigPlac.pdf")
