#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: Theo Portlock
For diversity analysis
'''
import functions as f
import pandas as pd
import seaborn as sns

# STOPPED HERE 24072023

# Baseline
base = pd.read_csv("../data/Donorandpatientbaselinecharacteristics.tsv", sep='\t', index_col=0)
base.drop([
    'Date-sample',
    'Hospitalisation in last 12 months', 
    'Travel in last 12 months'

# Meta
meta = pd.read_csv("../data/meta.csv", index_col=0)
meta['ARM'] = meta.Type + '_d' + meta['Days after treatment'].astype(int).astype(str)
var = meta.groupby('ARM').nunique().index
colours = pd.Series(sns.color_palette("Set2", len(var)).as_hex(), index=var)
colours.to_csv('../results/colours.csv')
meta.to_csv('../results/meta.tsv', sep='\t')

# Gut species
guttaxo = pd.read_csv(f"../data/gutTaxo.csv", sep='\t', index_col=0)
subject='gutmsp'
df = pd.read_csv(f"../data/{subject}.csv", index_col=0).T
df = df.join(meta['ARM'], how='inner').reset_index().set_index(['index', 'ARM'])
df = df.loc[df.sum(axis=1) != 0, df.sum(axis=0) != 0]
df.to_csv(f'../results/{subject}.csv')
df = f.norm(df).droplevel(1)
df = df.T.join(guttaxo['species']).groupby('species').sum().T
df.to_csv(f'../results/{subject}.tsv', sep='\t')

# Oral species
subject='oralmsp'
taxo = pd.read_csv(f"../data/oraltaxo.csv", index_col=0)
df = pd.read_csv(f"../data/{subject}.csv", index_col=0).T
df = df.join(meta['ARM'], how='inner').reset_index().set_index(['index', 'ARM'])
df = df.loc[df.sum(axis=1) != 0, df.sum(axis=0) != 0]
df.to_csv(f'../results/{subject}.csv')
df = f.norm(df).droplevel(1)
df = df.T.join(taxo['species']).groupby('species').sum().T
df.to_csv(f'../results/{subject}.tsv', sep='\t')

# need to update after this
# Pathways species
subject='pathway_community'
df = pd.read_csv(f"../results/{subject}.tsv", sep='\t', index_col=[0,1])
diversityanalysis(df, subject)
plotting(subject)

# Metabolomics
subject='metab'
df = pd.read_csv(f"../results/{subject}.csv", index_col=[0,1])
diversityanalysis(df, subject)
plotting(subject)
