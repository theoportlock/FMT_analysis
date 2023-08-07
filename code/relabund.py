#!/usr/bin/env python
import functions
import pandas as pd
import matplotlib.pyplot as plt

meta = pd.read_csv("../../data/meta.csv").set_index('Sample ID')
gutmsp = pd.read_csv("../../data/gutmsp.csv", index_col=0).T
oralmsp = pd.read_csv("../../data/oralmsp.csv", index_col=0).T
guttaxo = pd.read_csv("../../data/guttaxo.csv", index_col=0)
oraltaxo = pd.read_csv("../../data/oraltaxo.csv", index_col=0)

x = 'Type'
taxoType='species'
variables = ['Days after treatment', 'Type', 'Patient_Id']

taxoOralmsp = oralmsp.T.join(oraltaxo[taxoType]).groupby(taxoType).sum().T
metaTaxoOralmsp = taxoOralmsp.join(meta[variables], how='inner').set_index(variables)
functions.abund(metaTaxoOralmsp.sort_index())
plt.show()

taxoGutmsp = gutmsp.T.join(guttaxo[taxoType]).groupby(taxoType).sum().T
metaTaxoGutmsp = taxoGutmsp.join(meta[variables], how='inner').set_index(variables)
functions.abund(metaTaxoGutmsp.sort_index())
plt.show()

cor, pval = functions.corr(metaTaxoGutmsp.add_prefix('gut '), metaTaxoOralmsp.add_prefix('oral '))
functions.clustermap(cor,pval)
