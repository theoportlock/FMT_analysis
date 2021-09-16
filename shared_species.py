#!/usr/bin/env python
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy.stats import mannwhitneyu

samples_metadata = pd.read_csv('newmergedmetadata.csv').drop_duplicates(subset='ID_x', keep='first').set_index('ID_x').sort_index().dropna(thresh=20, axis=1)
samples_sm = pd.read_csv("antismashNorm.csv", index_col=0).T

samples_metadata.loc[samples_metadata.Type == 'DONOR', 'Days after treatment'] = 0

gmsp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
omsp_samples = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)
gmsp_gtaxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
omsp_otaxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)

taxaType='species'
gtaxonomy_samples = gmsp_samples.join(gmsp_gtaxonomy[taxaType], how='inner').groupby(taxaType).sum() > 0
samples_gtaxonomy = gtaxonomy_samples.T

otaxonomy_samples = omsp_samples.join(omsp_otaxonomy[taxaType], how='inner').groupby(taxaType).sum() > 0
samples_otaxonomy = otaxonomy_samples.T

sharedind = pd.Index.intersection(samples_gtaxonomy.index, samples_otaxonomy.index)
sharedcol = pd.Index.intersection(samples_gtaxonomy.columns, samples_otaxonomy.columns)

sharedvals = samples_gtaxonomy.loc[sharedind, sharedcol] & samples_otaxonomy.loc[sharedind, sharedcol]
sharedvals_metadata = samples_otaxonomy.join(samples_metadata[['Days after treatment','Type']], how='inner')

sharedvals_metadata = samples_otaxonomy.join(samples_metadata[['Days after treatment','Type','Patient_ID']], how='inner')
plotdf = sharedvals_metadata.groupby(['Days after treatment', 'Type']).mean().reset_index()

# need tomelt = plotdf.melt(id_vars=['Days after treatment', 'Type']) do the melt

#melt = plotdf.melt(id_vars=['Days after treatment', 'Type'])
#fmelt =melt.set_index(['Days after treatment', 'Type', 'variable'])
#filterplotdf = plotdf.loc[[plotdf > 0.95]]
#filtmelt = fmelt.loc[fmelt.value > 0.8].reset_index()

#sharedvals_metadata = samples_otaxonomy.join(samples_metadata[['Days after treatment','Type','Patient_ID']], how='inner')
#sharedvals_metadata.loc[sharedvals_metadata.Type == 'FMT'].set_index(['Days after treatment', 'Patient_ID']).sum(axis=1).reset_index().groupby('Days after treatment').mean()
#sharedvals_metadata.loc[sharedvals_metadata.Type == 'PLACEBO'].set_index(['Days after treatment', 'Patient_ID']).sum(axis=1).reset_index().groupby('Days after treatment').mean()

pl = sharedvals_metadata\
        .loc[sharedvals_metadata.Type == 'PLACEBO']\
        .drop('Type', axis=1)\
        .sum(axis=1)
        .set_index('Days after treatment')\
        .sum(axis=1)\
        .reset_index()\
        .groupby('Days after treatment')\
        .agg([np.mean,np.std])\
        .droplevel(0, axis=1)

samples_gutsum = pd.DataFrame(gtaxonomy_samples.agg(np.count_nonzero, axis=0), columns=['Species richness']).join(samples_metadata)
samples_oralsum = pd.DataFrame(otaxonomy_samples.copy().agg(np.count_nonzero, axis=0), columns=['Species richness']).join(samples_metadata)

timepoint = 0
orsum = samples_oralsum\
        .loc[samples_oralsum.Type == 'FMT']\
        .groupby('Days after treatment')\
        .agg([np.mean,np.std])['Species richness']
gusum = samples_gutsum\
        .loc[samples_gutsum.Type == 'FMT']\
        .groupby('Days after treatment')\
        .agg([np.mean,np.std])['Species richness']
venn2(subsets=(gusum['mean'][timepoint], orsum['mean'][timepoint], pl['mean'][timepoint]), set_labels=('gut', 'oral'))

samples_otaxonomy.join(samples_metadata[['Days after treatment','Type']], how='inner').query("Type == 'FMT'").set_index('Days after treatment').agg(np.count_nonzero, axis=1).reset_index().groupby('Days after treatment').agg([np.mean, np.std]).droplevel(0, axis=1)
samples_gtaxonomy.join(samples_metadata[['Days after treatment','Type']], how='inner').query("Type == 'FMT'").set_index('Days after treatment').agg(np.count_nonzero, axis=1).reset_index().groupby('Days after treatment').agg([np.mean, np.std]).droplevel(0, axis=1)

mannwhitneyu(comparethese.query("Type == 'FMT' and `Days after treatment` == 0"), comparethese.query("Type == 'FMT' and `Days after treatment` == 7"))
