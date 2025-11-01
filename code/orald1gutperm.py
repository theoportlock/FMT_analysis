#!/usr/bin/env python
from scipy.spatial import distance
import functions as f
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

meta = pd.read_csv('../data/coremeta.csv', index_col=0)
meta = pd.read_csv('../data/meta.csv', index_col=0)
gmsp = pd.read_csv("../data/gutmsp.csv", index_col=0)
omsp = pd.read_csv("../data/oralmsp.csv", index_col=0)
gtaxo = pd.read_csv("../data/guttaxo.csv", index_col=0)
otaxo = pd.read_csv("../data/oraltaxo.csv", index_col=0)

taxaType='genus'

gtaxoMsp = gmsp.join(gtaxo[taxaType], how='inner').groupby(taxaType).sum()
otaxoMsp = omsp.join(otaxo[taxaType], how='inner').groupby(taxaType).sum()
gtaxoMsp = f.norm(gtaxoMsp).dropna()
otaxoMsp = f.norm(otaxoMsp).dropna()

sharedind = pd.Index.intersection(gtaxoMsp.index, otaxoMsp.index)
sharedcol = pd.Index.intersection(gtaxoMsp.columns, otaxoMsp.columns)

gtaxoMsp = gtaxoMsp[sharedcol]
otaxoMsp = otaxoMsp[sharedcol]
gtaxoMsp['tmp1'] = ''
otaxoMsp['tmp2'] = ''
fgtaxoMsp = gtaxoMsp.join(otaxoMsp['tmp2'],how='outer').fillna(0.0)
fotaxoMsp = otaxoMsp.join(gtaxoMsp['tmp1'],how='outer').fillna(0.0)

gtaxoMsp = fgtaxoMsp.drop(['tmp1','tmp2'], axis=1)
otaxoMsp = fotaxoMsp.drop(['tmp1','tmp2'], axis=1)
nmeta = meta.loc[otaxoMsp.columns]

brays=pd.Series()
indiv = sharedcol[0]
for indiv in sharedcol:
    try:
        oindiv = meta.loc[meta.id == meta.loc[indiv].id].query('`Days after treatment` == 0').index
        brays[indiv] = distance.braycurtis(gtaxoMsp[indiv], otaxoMsp[oindiv].iloc[:,0])
    except:
        pass


out = brays.to_frame('Bray Curtis').join(meta[['Type','Days after treatment']], how='inner').set_index(['Type','Days after treatment'])
f.setupplot()
f.box(df=out.reset_index(),
      x='Days after treatment',
       y='Bray Curtis',
       hue='Type'
       )
plt.savefig('../results/betaspeciesoralbaseline.svg')

f.change(out.xs('FMT')).T
