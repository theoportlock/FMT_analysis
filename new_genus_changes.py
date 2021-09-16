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
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set(font_scale=0.8)

samples_metadata = pd.read_csv("newmergedmetadata.csv").set_index('ID_x')
#samples_metadata = pd.read_csv('metadata.csv', index_col=0)

taxaType='genus'
gmsp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
gmsp_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
samples_gtaxonomy = gmsp_samples.join(gmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
gmspdf = samples_gtaxonomy.add_prefix('gut ')

omsp_samples = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)
omsp_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)
samples_otaxonomy = omsp_samples.join(omsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
omspdf = samples_otaxonomy.add_prefix('oral ')

variable=['Days after treatment', 'Type']
df = omspdf.join(gmspdf, how='inner')
df = df.join(samples_metadata[variable], how='inner')
#meandf = df.loc[df.Type == 'FMT'].groupby('Days after treatment').mean()
#baseline = df.loc[df['Days after treatment'] == 0]
#baseline = meandf.loc[0]
#notbaseline  = df.loc[df['Days after treatment'] != 0]
average = df.groupby('Days after treatment').mean()
baseline = df.loc[(df.Type == 'FMT') & (df['Days after treatment'] == 0)].drop("Type", axis =1)
difference = average.sub(baseline.mean())
#baseline = meandf.loc[0]
notbaseline  = df.loc[(df.Type == 'FMT') & (df['Days after treatment'] != 0)].drop("Type", axis =1)
#subdf = meandf.sub(baseline).T.drop(0, axis=1) * 1000000
#subdf = meandf.sub(baseline).T.drop(0, axis=1)
#subdf.to_csv('subdf.csv')

def mandf(df):
    #base = baseline.loc[(baseline.sum(axis=1)!=0),(baseline.sum(axis=0)!=0)].copy()
    #day = df.loc[(df.sum(axis=1)!=0),(df.sum(axis=0)!=0)].copy()
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
#stats = notbaseline.groupby('Days after treatment').apply(mandf)
                
significantMatrix = stats.loc[:, (stats < 0.05).any(axis=0)].drop('Days after treatment', axis=1)
plotMatrix = difference[significantMatrix.columns]
plotMatrix = plotMatrix.T
significantMatrix = significantMatrix.T.iloc[:,1:]
significantMatrix = significantMatrix.loc[plotMatrix.index]

