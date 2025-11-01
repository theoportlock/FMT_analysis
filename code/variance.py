#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.spatial import distance
from skbio.stats.distance import permanova
import argparse
import numpy as np
import pandas as pd
import skbio
import functions as f

def PERMANOVA(df, pval=True, full=False):
    Ar_dist = distance.squareform(distance.pdist(df, metric="braycurtis"))
    DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
    result = permanova(DM_dist, df.index)
    if full:
        return result
    if pval:
        return result['p-value']
    else:
        return result['test statistic']

meta = pd.read_csv("../data/cleanMeta.csv").set_index('Sample ID')

oralmsp = pd.read_csv("../data/oralmsp.csv", index_col=0)
oraltaxo = pd.read_csv("../data/oraltaxo.csv", index_col=0)
oralmsp = f.norm(oralmsp.join(oraltaxo['species']).groupby('species').sum().T)
oralcard = pd.read_csv("../data/oralcard.csv", index_col=0).T
oralcazy = pd.read_csv("../data/oralcazy.csv", index_col=0).T
oralkegg = pd.read_csv("../data/oralkegg.csv", index_col=0).T
oralpfam = pd.read_csv("../data/oralpfam.csv", index_col=0).T

gutmsp = pd.read_csv("../data/gutmsp.csv", index_col=0)
guttaxo = pd.read_csv("../data/guttaxo.csv", index_col=0)
gutmsp = f.norm(gutmsp.join(guttaxo['species']).groupby('species').sum().T)
gutcard = pd.read_csv("../data/gutcard.csv", index_col=0).T
gutcazy = pd.read_csv("../data/gutcazy.csv", index_col=0).T
gutkegg = pd.read_csv("../data/gutkegg.csv", index_col=0).T
gutpfam = pd.read_csv("../data/gutpfam.csv", index_col=0).T

cytokines = pd.read_csv("../data/cytokines.csv", index_col=0)

plasmaBA = pd.read_csv("../data/plasmaBileAcid.tsv", sep='\t', index_col=0)
plasmaTrp = pd.read_csv("../data/plasmaTryptophan.tsv", sep='\t', index_col=0)
urineBA = pd.read_csv("../data/urineBileAcid.tsv", sep='\t', index_col=0)
urineNMR = pd.read_csv("../data/urineNmr.tsv", sep='\t', index_col=0)
urineTrp = pd.read_csv("../data/urineTryptophan.tsv", sep='\t', index_col=0)
stoolBA = pd.read_csv("../data/stoolBileAcid.tsv", sep='\t', index_col=0)
stoolNMR = pd.read_csv("../data/stoolNmr.tsv", sep='\t', index_col=0)
stoolTrp = pd.read_csv("../data/stoolTryptophan.tsv", sep='\t', index_col=0)

# filter meta
meta = meta.loc[meta.Type != 'DONOR']

datasets = {'oralmsp': oralmsp,
'oralcard ': oralcard ,
'oralcazy ': oralcazy ,
'oralkegg ': oralkegg ,
'oralpfam ': oralpfam ,
'gutmsp ': gutmsp ,
'gutcard ': gutcard ,
'gutcazy ': gutcazy ,
'gutkegg ': gutkegg ,
'gutpfam ': gutpfam ,
'cytokines ': cytokines ,
'plasmaBA ': plasmaBA ,
'plasmaTrp ': plasmaTrp ,
'urineBA ': urineBA ,
'urineNMR ': urineNMR ,
'urineTrp ': urineTrp ,
'stoolBA ': stoolBA ,
'stoolNMR ': stoolNMR ,
'stoolTrp': stoolTrp}

output = pd.DataFrame(index=meta['Days after treatment'].sort_values().unique(), 
                        columns=datasets.keys())
for name, dataset in datasets.items():
    joined = dataset.join(meta[['Days after treatment', 'Type']], how='inner')
    #time = joined['Days after treatment'].sort_values().unique()[0]
    for time in joined['Days after treatment'].sort_values().unique():
        tomeasure = joined.loc[joined['Days after treatment'] == time].drop('Days after treatment', axis=1).set_index('Type')
        output.loc[time, name] = PERMANOVA(tomeasure)

power = -output.astype(float).apply(np.log).T
sig = -np.log(0.05)
f.setupplot()
heatmap(power, power.gt(sig))
