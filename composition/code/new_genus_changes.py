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
from statsmodels.stats.multitest import fdrcorrection
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set(font_scale=0.8)

samples_metadata = pd.read_csv("../../data/newnewmetadata.csv").set_index('Sample ID')
samples_gutmsp = pd.read_csv("../../data/gutmsp.csv", index_col=0)
samples_oralmsp = pd.read_csv("../../data/oralmsp.csv", index_col=0)

taxaType='genus'
gmsp_taxonomy = pd.read_csv("../../data/guttaxo.csv", index_col=0)
samples_gtaxonomy = samples_gutmsp.join(gmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
gmspdf = samples_gtaxonomy.add_prefix('gut ')
omsp_taxonomy = pd.read_csv("../../data/oraltaxo.csv", index_col=0)
samples_otaxonomy = samples_oralmsp.join(omsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
omspdf = samples_otaxonomy.add_prefix('oral ')

variable=['Days after treatment', 'Type']
df = omspdf.join(gmspdf, how='inner')
df = df.join(samples_metadata[variable], how='inner')
average = df.groupby('Days after treatment').mean()
baseline = df.loc[(df.Type == 'FMT') & (df['Days after treatment'] == 0)].drop("Type", axis =1)
#difference = average.sub(baseline.mean()).div(baseline.mean())
difference = average.div(baseline.mean())
notbaseline  = df.loc[(df.Type == 'FMT') & (df['Days after treatment'] != 0)].drop("Type", axis =1)

def mandf(df):
    base = baseline
    day = df
    result = pd.Series()
    for i in base.columns:
        try:
            result[i] = mannwhitneyu(base[i],day[i])[1]
        except:
            result[i] = np.nan
    return result

stats = notbaseline.groupby('Days after treatment').apply(mandf)
significantMatrix = stats.loc[:, (stats < 0.05).any(axis=0)].drop('Days after treatment', axis=1)

significantMatrix = pd.DataFrame(
    fdrcorrection(stats.values.flatten())[0].reshape(stats.shape),
    index = stats.index,
    columns = stats.columns)

plotMatrix = difference[significantMatrix.columns]
plotMatrix = plotMatrix.T
significantMatrix = significantMatrix.T
significantMatrix = significantMatrix.loc[plotMatrix.index]

