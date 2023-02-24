#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: Theo Portlock
For diversity analysis
'''
import functions as f
import pandas as pd
import seaborn as sns

meta = pd.read_csv("../data/meta.csv", index_col=0)
meta['ARM'] = meta.Type + '_d' + meta['Days after treatment'].astype(int).astype(str)
var = meta.groupby('ARM').nunique().index
colours = pd.Series(sns.color_palette("Set2", len(var)).as_hex(), index=var)
colours.to_csv('../results/colours.csv')

subject='gutmsp'
df = pd.read_csv(f"../data/{subject}.csv", index_col=0).T
df = df.join(meta['ARM'], how='inner').reset_index().set_index(['index', 'ARM'])

subject='oralmsp'
df = pd.read_csv(f"../data/{subject}.csv", index_col=0).T
df =df.join(meta.ARM, how='inner').set_index('ARM')
diversityanalysis(df, subject)
plotting(subject)

subject='pathway_community'
df = pd.read_csv(f"../results/{subject}.tsv", sep='\t', index_col=[0,1])
diversityanalysis(df, subject)
plotting(subject)

subject='metab'
df = pd.read_csv(f"../results/{subject}.csv", index_col=[0,1])
diversityanalysis(df, subject)
plotting(subject)
