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
meandf = df.loc[df.Type == 'PLACEBO'].groupby('Days after treatment').mean()
baseline = meandf.loc[0]
notbaseline  = df.loc[(df.Type == 'PLACEBO') & (df['Days after treatment'] != 0)].drop("Type", axis =1)
difference = meandf.sub(baseline.mean())
subdf = meandf.sub(baseline).T.drop(0, axis=1) * 1000000
#subdf.to_csv('subdf.csv')
