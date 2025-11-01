#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Plotting an annotated seaborn correlation heatmap between MSP data and metadata
Theo Portlock
'''

import seaborn as sns
import pandas as pd
import functions as f
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

meta = pd.read_csv('../data/cleanMeta.csv').set_index('Sample ID').fillna(0)._get_numeric_data()
meta.columns = meta.columns.str.replace(' ','_')
meta = meta['MELDPLUS'].to_frame()
othermeta = pd.read_csv('../data/plasmaBileAcid.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('plasma ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/plasmaTryptophan.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('plasma ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/stoolBileAcid.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('stool ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/stoolNmr.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('stool ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/stoolTryptophan.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('stool ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/urineBileAcid.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('urine ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/urineNmr.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('urine ')
meta = meta.join(othermeta)
othermeta = pd.read_csv('../data/urineTryptophan.tsv', sep='\t').set_index('Sample ID')
othermeta = othermeta.add_prefix('urine ')
meta = meta.join(othermeta)

meta = meta.loc[meta.MELDPLUS != 0]
meta.dropna(inplace=True)

taxo = pd.read_csv("../data/guttaxo.csv", index_col=0)
msp = pd.read_csv("../data/gutmsp.csv", index_col=0)
taxaType='species'
mspTaxo = msp.join(taxo[taxaType], how='inner').groupby(taxaType).sum().T
mspTaxo = mspTaxo.add_prefix('gut ')
meta = meta.join(mspTaxo, how='inner').astype('float')

taxo = pd.read_csv("../data/oraltaxo.csv", index_col=0)
msp = pd.read_csv("../data/oralmsp.csv", index_col=0)
taxaType='species'
mspTaxo = msp.join(taxo[taxaType], how='inner').groupby(taxaType).sum().T
mspTaxo = mspTaxo.add_prefix('oral ')
meta = meta.join(mspTaxo, how='inner').astype('float')

#kegg = pd.read_csv("../data/gutkegg.csv", index_col=0).T
#meta = meta.join(kegg).dropna()

meta = meta.loc[:,f.richness(meta.T).gt(10).Richness]
cor = meta.corr()
edges = f.to_edges(cor, thresh=0.4)
edges.to_csv('../results/edges.csv')

nodes, pval = f.corr(meta['MELDPLUS'].to_frame(), meta.drop(['MELDPLUS'], axis=1))
nodes.T.sort_values('MELDPLUS').to_csv('../results/nodes.csv')

