#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Compositional changes in the FMT dataset
Theo Portlock
'''
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import functions as f

meta = pd.read_csv("../data/newnewmetadata.csv").set_index('Sample ID')
gutmsp = pd.read_csv("../data/gutmsp.csv", index_col=0)
oralmsp = pd.read_csv("../data/oralmsp.csv", index_col=0)

taxaType='species'
gmsp_taxonomy = pd.read_csv("../data/gutTaxo.csv", index_col=0, sep='\t')
gutmsp = f.norm(gutmsp.T)
gtaxonomy = f.taxofunc(gutmsp, gmsp_taxonomy, short=True)
gtaxonomy = gtaxonomy.loc[~gtaxonomy.index.isna()].T
#gtaxonomy = gutmsp.join(gmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
gmspdf = gtaxonomy.add_prefix('gut ')

'''
omsp_taxonomy = pd.read_csv("../data/oraltaxo.csv", index_col=0)
otaxonomy = oralmsp.join(omsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
omspdf = otaxonomy.add_prefix('oral ')
omspdf = f.minmaxscale(omspdf)
'''

variable=['Days after treatment', 'Type']
df = omspdf.join(gmspdf, how='inner')
df = gmspdf.copy()
df = df.join(meta[variable], how='inner').set_index(variable)
df = f.mult(f.filter(df, min_unique=4))
fc = f.lfc(df).T
sig = f.sig(df).T

subject='species'
fc.to_csv(f'../results/{subject}FCabundance.csv')
sig.to_csv(f'../results/{subject}MWWabundance.csv')
fcthresh=0
val=0.0005

ffc = fc.loc[:, fc.abs().gt(fcthresh).any(axis=0) & sig.lt(val).any(axis=0)]
fsig = sig.loc[:, fc.abs().gt(fcthresh).any(axis=0) & sig.lt(val).any(axis=0)]
f.clustermap(ffc, fsig.lt(val), figsize=(3,7))
plt.savefig(f'../results/{subject}ClusterFC_HvD.svg')

