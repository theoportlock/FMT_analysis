#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Plotting an annotated seaborn correlation heatmap between MSP data and metadata
Theo Portlock
'''
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import functions as f

f.setupplot()
samples_metadata = pd.read_csv('../data/cleanMeta.csv').set_index('Sample ID').fillna(0)._get_numeric_data()
ometa = []
sampleslist=[
    '../data/plasmaBileAcid.tsv',
    '../data/plasmaTryptophan.tsv',
    '../data/stoolBileAcid.tsv', 
    '../data/stoolNmr.tsv',
    '../data/stoolTryptophan.tsv',
    '../data/urineBileAcid.tsv', 
    '../data/urineNmr.tsv',
    '../data/urineTryptophan.tsv']

othermeta = pd.read_csv('../data/plasmaBileAcid.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/plasmaTryptophan.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/stoolBileAcid.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/stoolNmr.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/stoolTryptophan.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/urineBileAcid.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/urineNmr.tsv', sep= '\t').set_index('Sample ID')
othermeta = pd.read_csv('../data/urineTryptophan.tsv', sep= '\t').set_index('Sample ID')

#samples_metadata.drop("Creatinine", axis=1, inplace=True)
#samples_metadata = samples_metadata.join(othermeta)
samples_metadata = samples_metadata.loc[samples_metadata.MELDNA != 0]
samples_metadata.dropna(inplace=True)
msp_taxonomy = pd.read_csv("../data/guttaxo.csv", index_col=0)
msp_samples = pd.read_csv("../data/gutmsp.csv", index_col=0)
msp_taxonomy = pd.read_csv("../data/oraltaxo.csv", index_col=0)
msp_samples = pd.read_csv("../data/oralmsp.csv", index_col=0)

# Join taxo information
taxaType='genus'
taxonomy_samples = msp_samples.join(msp_taxonomy[taxaType], how='inner').groupby(taxaType).sum()
samples_taxonomy = taxonomy_samples.T
samples_taxonomyMetadata = samples_taxonomy.join(samples_metadata, how='inner').astype('float')

# Calculate and format correlations
thresh = 0.005
cor, pval = f.corr(samples_taxonomy, samples_metadata, min_unique=10)
cor = cor.loc[pval.lt(thresh).sum(axis=1) != 0, pval.lt(thresh).sum(axis=0) != 0]
pval = pval.loc[pval.lt(thresh).sum(axis=1) != 0, pval.lt(thresh).sum(axis=0) != 0]

# Plot
f.clustermap(cor, pval.lt(thresh), figsize=(4,4))
f.clustermap(cor, pval.lt(thresh), figsize=(3,1.5))
plt.savefig("../results/gutmetaclusterclustermap.svg")
plt.savefig("../results/oralmetaclusterclustermap.svg")
plt.show()
