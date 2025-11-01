#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: Theo Portlock
For diversity analysis
'''
import pandas as pd

guttaxo = pd.read_csv(f"data/gutTaxo.csv", sep='\t', index_col=0)
species = (guttaxo.species + ':' + guttaxo.index).to_frame('species')

subject='gutmsp'
df = pd.read_csv(f"data/{subject}.csv", index_col=0).T
df = df.loc[df.sum(axis=1) != 0, df.sum(axis=0) != 0]
df = df.T.div(df.sum(axis=1), axis=1).T

df = df.T.join(species).dropna().set_index('species').T

df.to_csv(f'results/species.tsv', sep='\t')
