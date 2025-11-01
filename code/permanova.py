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

meta = f.load('meta')
gutmsp = pd.read_csv("../data/gutmsp.csv", index_col=0)

taxaType='species'
gmsp_taxonomy = pd.read_csv("../data/gutTaxo.csv", index_col=0, sep='\t')
gutmsp = f.norm(gutmsp.T)
gtaxonomy = gutmsp.T.join(gmsp_taxonomy['species'], how='inner').set_index('species').T

meta = meta.loc[meta['Timepoint'] == 0.0]
meta['Type'] = meta.Type.map({'FMT':'Participant', 'PLACEBO':'Participant'}).fillna('Donor')

strat = f.stratify(gtaxonomy, meta, 'Type')

pcoa = f.calculate('pcoa', strat)
 
f.setupplot(figsize=(3,3))
f.spindle(pcoa)
plt.xlabel('PCoA1 (18.5%)')
plt.ylabel('PCoA1 (13.3%)')
f.savefig('spindle')
f.PERMANOVA(strat)

meta = f.load('meta')
meta = meta.loc[meta.Type != 'DONOR']

strat = gtaxonomy.join(meta[['id','Timepoint','Type']], how='inner').set_index(['id','Timepoint','Type'])

from itertools import permutations
from scipy.spatial import distance
import skbio

combs = list(permutations(strat.index.unique(), 2))
outdf = pd.DataFrame(index=pd.MultiIndex.from_tuples(combs), columns=['beta','pval'])

df = strat.copy()
Ar_dist = distance.squareform(distance.pdist(df, metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
betas = pd.DataFrame(DM_dist.data, index=df.index, columns=df.index)
baselinebetas = betas.xs(0.0, level=1, axis=1).droplevel(1, axis=1)
stacked = baselinebetas.stack()

out = []
for i in stacked.index:
    if i[0] == i[3]:
        out.append(i)

fstacked = stacked.loc[out].to_frame('Beta diversity').droplevel(0).reset_index()
f.setupplot(figsize=(3,3))
f.box(df = fstacked, x='Timepoint', y='Beta diversity', hue='Type')
f.savefig('betadiv')
fstacked.index = fstacked.Type +'_' + fstacked.Timepoint.astype(int).astype(str)
fstacked = fstacked['Beta diversity'].to_frame()
ch = f.change(fstacked)
t1 = ch['PLACEBO_7vsFMT_7']
t2 = ch['PLACEBO_30vsFMT_30']
t3 = ch['PLACEBO_90vsFMT_90']
