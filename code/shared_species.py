#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns

meta = pd.read_csv('../data/coremeta.csv', index_col=0)
gmsp = pd.read_csv("../data/gutmsp.csv", index_col=0)
omsp = pd.read_csv("../data/oralmsp.csv", index_col=0)
gtaxo = pd.read_csv("../data/guttaxo.csv", index_col=0)
otaxo = pd.read_csv("../data/oraltaxo.csv", index_col=0)

taxaType='genus'
v = ['Days after treatment','Type','Patient_Id']

gtaxoMsp = gmsp.join(gtaxo[taxaType], how='inner').groupby(taxaType).sum() > 0
otaxoMsp = omsp.join(otaxo[taxaType], how='inner').groupby(taxaType).sum() > 0

sharedind = pd.Index.intersection(gtaxoMsp.index, otaxoMsp.index)
sharedcol = pd.Index.intersection(gtaxoMsp.columns, otaxoMsp.columns)

shared = gtaxoMsp.loc[sharedind, sharedcol] & otaxoMsp.loc[sharedind, sharedcol]

metaShared = shared.T.join(meta[v], how='inner').set_index(v)

probaMetaShared = metaShared.groupby(level=[0,1]).mean().sort_index(level=1)
#absMetaShared = metaShared.groupby(level=[0,1]).mean().ge(0.9).sort_index(level=1)

probaMetaShared = probaMetaShared.T[probaMetaShared.sum() != 0]
probaMetaShared = probaMetaShared.loc[~probaMetaShared.index.str.contains('unclassified')]

f.clustermap(probaMetaShared)
plt.savefig('../results/sharedspecies.pdf')
plt.show()
