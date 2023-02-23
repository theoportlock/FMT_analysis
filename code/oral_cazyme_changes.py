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
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set(font_scale=0.8)

samples_metadata = pd.read_csv("../../data/newnewmetadata.csv").set_index('Sample ID')
data_samples = pd.read_csv("../../data/oralcazy.csv", index_col=0)
samples_data=data_samples.T

# Join data information
samples_dataMetadata = samples_data.join(samples_metadata['Days after treatment'], how='inner')

# Calculate and format changes
averageDays = samples_dataMetadata.groupby('Days after treatment').mean()
baseline = samples_dataMetadata.loc[samples_dataMetadata['Days after treatment'] == 0]
#difference = averageDays.sub(baseline.mean())
difference = averageDays.div(baseline.mean()).apply(np.log)
difference = difference.loc[(difference.sum(axis=1) != 0), (difference.sum(axis=0) != 0)]
difference = difference.loc[(difference.sum(axis=1) != 1), (difference.sum(axis=0) != 1)]
difference.replace([np.inf, -np.inf], np.nan, inplace=True)
difference.dropna(axis=1, inplace=True)

def mandf(df):
    base = baseline.loc[(baseline.sum(axis=1)!=0),(baseline.sum(axis=0)!=0)].copy()
    day = df.loc[(df.sum(axis=1)!=0),(df.sum(axis=0)!=0)].copy()
    result = pd.Series()
    for i in base.columns:
        try:
            result[i] = mannwhitneyu(base[i],day[i])[1]
        except:
            result[i] = np.nan
    return result

stats = samples_dataMetadata.groupby('Days after treatment').apply(mandf)
                
#slicedUncorrectedPValues = stats.loc[:, (stats < 0.05).any(axis=0)]
significantMatrix = stats.loc[:, (stats < 0.05).any(axis=0)]
#significantMatrix = pd.DataFrame(
#    fdrcorrection(slicedUncorrectedPValues.values.flatten())[0].reshape(slicedUncorrectedPValues.shape),
#    index = slicedUncorrectedPValues.index,
#    columns = slicedUncorrectedPValues.columns)
plotMatrix = difference[significantMatrix.columns]
plotMatrix = plotMatrix.T.iloc[:,1:]
significantMatrix = significantMatrix.T.iloc[:,1:]
significantMatrix = significantMatrix.loc[plotMatrix.index]

#plotMatrix = values.loc[significantMatrix.index]

# Plot
g = sns.clustermap(
    data=plotMatrix,
    col_cluster=False,
    row_cluster=False,
    cmap="coolwarm",
#    vmax=50, 
#    vmin=-50, 
    yticklabels=True,
    xticklabels=True)

for tick in g.ax_heatmap.get_yticklabels():
    tick.set_rotation(0)

annot = pd.DataFrame(index=plotMatrix.index, columns=plotMatrix.columns)
annot[(significantMatrix < 0.0005) & (plotMatrix > 0)] = '+++'
annot[(significantMatrix < 0.005) & (plotMatrix > 0)] = '++'
annot[(significantMatrix < 0.05) & (plotMatrix > 0)] = '+'
annot[(significantMatrix < 0.0005) & (plotMatrix < 0)] = '---'
annot[(significantMatrix < 0.005) & (plotMatrix < 0)] = '--'
annot[(significantMatrix < 0.05) & (plotMatrix < 0)] = '-'
annot[significantMatrix >= 0.05] = ''

for i, ix in enumerate(plotMatrix.index):
    for j, jx in enumerate(plotMatrix.columns):
        text = g.ax_heatmap.text(
            j + 0.5,
            i + 0.5,
            annot.values[i,j],
            ha="center",
            va="center",
            color="black",
        )
        text.set_fontsize(10)

plt.savefig("../results/oral_cazyme_changes.pdf")
#plt.show()
