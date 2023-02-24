#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: Theo Portlock
For diversity analysis
'''
import functions as f
import itertools
import matplotlib.pyplot as plt
import pandas as pd

def diversityanalysis(df, subject):
    diversity = pd.concat(
            [f.evenness(df),
             f.richness(df),
             f.shannon(df)],
            axis=1).sort_index().sort_index(ascending=False)
    diversity.to_csv(f'../results/{subject}Diversity.csv')
    diversity = diversity.droplevel(0)
    fc = f.lfc(diversity)
    fc.columns = fc.columns.str.join('/')
    fc.T.to_csv(f'../results/{subject}FCdiversity.csv')
    sig = f.sig(diversity)
    sig.columns = sig.columns.str.join('/')
    sig.T.to_csv(f'../results/{subject}MWWdiversity.csv')
    ndf = df.droplevel(0)
    pcoa = f.PCOA(ndf.sort_index())
    pcoa.to_csv(f'../results/{subject}PCOA.csv')
    permanova = f.PERMANOVA(ndf.reset_index(drop=True), ndf.index)
    with open(f'../results/{subject}permanova.txt','w') as of: of.write(permanova.to_string())
    # pairwise
    comb = list(itertools.combinations(ndf.index.unique(), 2))
    out= pd.DataFrame()
    i = comb[0]
    for i in comb: 
        out[i] = f.PERMANOVA(ndf.loc[list(i)].reset_index(drop=True), ndf.loc[list(i)].index)
    out.to_csv(f'../results/{subject}pairwiseANOVA.csv')

def plotting(subject):
    f.setupplot()
    colours = pd.read_csv('../results/colours.csv', index_col=0)['0']
    diversity = pd.read_csv(f'../results/{subject}Diversity.csv', index_col=[0,1]).droplevel(0)
    metric = diversity.columns[0]
    stats = pd.read_csv(f'../results/{subject}MWWdiversity.csv', index_col=0)
    stats.index=stats.index.str.split('/')
    metric = diversity.columns[2]
    for metric in diversity.columns:
        fig, ax = plt.subplots(figsize=(1.5,1.5))
        kwargs = {
            'data':diversity,
            'x':diversity.index,
            'y':metric,
            'palette':colours.to_dict(),
            'dodge':False,
            'ax':ax,
            'stats':stats[metric].lt(0.05),
            }
        f.box(**kwargs)
        plt.savefig(f'../results/{subject}box{metric.split()[0]}.svg')
    pcoa = pd.read_csv(f'../results/{subject}PCOA.csv', index_col=0)
    fig, ax = plt.subplots(figsize=(3,2.818))
    f.spindleplot(pcoa, ax=ax, palette=colours)
    plt.xlabel('PCoA1')
    plt.ylabel('PCoA2')
    plt.tight_layout()
    plt.savefig(f'../results/{subject}PCOA.svg')


meta = pd.read_csv("../data/meta.csv", index_col=0)
meta['ARM'] = meta.Type + '_d' + meta['Days after treatment'].astype(int).astype(str)

subject='gutmsp'
df = pd.read_csv(f"../data/{subject}.csv", index_col=0).T
df = df.join(meta['ARM'], how='inner').reset_index().set_index(['index', 'ARM'])
diversityanalysis(df, subject)
plotting(subject)

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
