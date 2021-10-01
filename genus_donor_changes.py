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
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import mannwhitneyu
import itertools
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set(font_scale=0.8)

variables = ['Days after treatment', 'Type']
taxaType='genus'
samples_metadata = pd.read_csv('newnewmetadata.csv', index_col=0).set_index('Sample ID').iloc[1:]

gmsp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
gmsp_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
samples_gtaxonomy = gmsp_samples.join(gmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
gmspdf = samples_gtaxonomy.add_prefix('gut ')

df = gmspdf.join(samples_metadata[variables], how='inner')
df.fillna(0,inplace=True)

#box_pairs = list(itertools.combinations([(a,b) for a in df['Days after treatment'].unique() for b in df["Type"].unique()],2))
box_pairs = [
        ((0.0, 'DONOR'), (0.0, 'PLACEBO')),
        ((0.0, 'DONOR'), (7.0, 'PLACEBO')),
        ((0.0, 'DONOR'), (30.0, 'PLACEBO')),
        ((0.0, 'DONOR'), (90.0, 'PLACEBO')),
        ((0.0, 'DONOR'), (0.0, 'FMT')),
        ((0.0, 'DONOR'), (7.0, 'FMT')),
        ((0.0, 'DONOR'), (30.0, 'FMT')),
        ((0.0, 'DONOR'), (90.0, 'FMT'))]


test_short_name = 'M.W.W'
pvalues = pd.DataFrame(index=joined.columns)
for pair in box_pairs:
    data1 = df.groupby(variables).get_group(pair[0]).drop(variables, axis=1)
    data2 = df.groupby(variables).get_group(pair[1]).drop(variables, axis=1)
    result = pd.Series()
    for i in pvalues.index:
        try:
            result[i] = mannwhitneyu(data1[i],data2[i])[1]
        except:
            result[i] = 1
    pvalues[pair] = result

values = pd.DataFrame(index=joined.columns)
for pair in box_pairs:
    data1 = df.groupby(variables).get_group(pair[0]).drop(variables, axis=1)
    data2 = df.groupby(variables).get_group(pair[1]).drop(variables, axis=1)
    result = pd.Series()
    for i in values.index:
        try:
            result[i] = np.mean(data2[i]) - np.mean(data1[i])
        except:
            result[i] = 1
    values[pair] = result

#sigpvals = pvalues < 0.05
significantMatrix = pvalues.loc[(pvalues < 0.05).any(axis=1), : ]

plotMatrix = values.loc[significantMatrix.index]
#plotMatrix = plotMatrix.T

# Plot
g = sns.clustermap(
    data=np.tanh(plotMatrix*50000),
    col_cluster=False,
    row_cluster=False,
    cmap="vlag",
#    vmin=-1e-8,
#    vmax=1e-6,
    center=0,
    yticklabels=True,
    xticklabels=True)

for tick in g.ax_heatmap.get_yticklabels():
    tick.set_rotation(0)

annot=pd.DataFrame(columns=plotMatrix.columns, index=plotMatrix.index)
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
        #text = g.text(
            j + 0.5,
            i + 0.5,
            annot.values[i,j],
            ha="center",
            va="center",
            color="black",
        )
        text.set_fontsize(8)

plt.savefig("results/sigchangegenusdonor.pdf")
#plt.show()
