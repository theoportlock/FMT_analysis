#!/usr/bin/env python import numpy as np
import seaborn as sns
import numpy as np
import pandas as pd
import itertools
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

sns.set(rc={'figure.figsize':(11.7,30)})
plt.figure(figsize=(15,8))
sns.set(font_scale=0.8)

samples_metadata = pd.read_csv("newnewmetadata.csv").set_index('Sample ID').iloc[1:]
gut_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
gut_msp_data = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
oral_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)
oral_msp_data = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)

taxa_type='family'
x='Days after treatment'
y='Enterobacteriaceae'
hue='Type'
variables = ['Days after treatment', 'Type']
taxaType='genus'

dtaxa_not_sum = gut_msp_data.join(gut_taxonomy[taxa_type], how='inner').set_index(taxa_type).T
otaxa_not_sum = oral_msp_data.join(oral_taxonomy[taxa_type], how='inner').set_index(taxa_type).T
otaxa = otaxa_not_sum.groupby(otaxa_not_sum.columns, axis=1).sum()
dtaxa = dtaxa_not_sum.groupby(dtaxa_not_sum.columns, axis=1).sum()

otaxa['Site'], dtaxa['Site'] = 'Oral', 'Gut'
joined = pd.concat([otaxa, dtaxa])
df = joined.join(samples_metadata[['Days after treatment', 'Type', 'id']], how='inner')

meandf = df.loc[df.Type == 'PLACEBO'].groupby('Days after treatment').mean()
baseline = meandf.loc[0]
notbaseline  = df.loc[(df.Type == 'PLACEBO') & (df['Days after treatment'] != 0)].drop("Type", axis =1)
difference = np.log(meandf) - np.log(baseline.mean())
subdf = meandf.sub(baseline).T.drop(0, axis=1)
#subdf.to_csv('subdf.csv')

#box_pairs = list(itertools.combinations([(a,b) for a in df['Days after treatment'].unique() for b in df["Type"].unique()],2))
box_pairs = [
        ((0.0, 'PLACEBO'), (7.0, 'PLACEBO')),
        ((0.0, 'PLACEBO'), (30.0, 'PLACEBO')),
        ((0.0, 'PLACEBO'), (90.0, 'PLACEBO')),
        ((0.0, 'FMT'), (7.0, 'FMT')),
        ((0.0, 'FMT'), (30.0, 'FMT')),
        ((0.0, 'FMT'), (90.0, 'FMT'))]


test_short_name = 'M.W.W'
pvalues = pd.DataFrame(index=joined.columns)
for pair in box_pairs:
    data1 = df.groupby(variables).get_group(pair[0]).drop(variables, axis=1)
    data2 = df.groupby(variables).get_group(pair[1]).drop(variables, axis=1)
    result = pd.Series()
    for i in pvalues.index:
        try:
            result[i] = mannwhitneyu(data1[i],data2[i])[1]
        except:
            result[i] = 1
    pvalues[pair] = result

values = pd.DataFrame(index=joined.columns)
for pair in box_pairs:
    data1 = df.groupby(variables).get_group(pair[0]).drop(variables, axis=1)
    data2 = df.groupby(variables).get_group(pair[1]).drop(variables, axis=1)
    result = pd.Series()
    for i in values.index:
        try:
            result[i] = np.mean(data2[i]) - np.mean(data1[i])
        except:
            result[i] = 1
    values[pair] = result

#significantMatrix = pvalues.loc[(pvalues < 0.05).any(axis=1), : ]
significantMatrix = pvalues

plotMatrix = values.loc[significantMatrix.index]

# Plot
g = sns.clustermap(
    data=np.tanh(plotMatrix*100000),
    col_cluster=False,
    row_cluster=False,
    cmap="vlag",
#    vmin=-1e-8,
#    vmax=1e-6,
    center=0,
    yticklabels=True,
    xticklabels=True)
g.fig.set_figwidth(8.27)
g.fig.set_figheight(11.7)
sns.set(font_scale=0.5)

for tick in g.ax_heatmap.get_yticklabels():
    tick.set_rotation(0)

annot=pd.DataFrame(columns=plotMatrix.columns, index=plotMatrix.index)
annot[(significantMatrix < 0.0005) & (plotMatrix > 0)] = '+++'
annot[(significantMatrix < 0.005) & (plotMatrix > 0)] = '++'
annot[(significantMatrix < 0.05) & (plotMatrix > 0)] = '+'
annot[(significantMatrix < 0.0005) & (plotMatrix < 0)] = '---'
annot[(significantMatrix < 0.005) & (plotMatrix < 0)] = '--'
annot[(significantMatrix < 0.05) & (plotMatrix < 0)] = '-'
annot[significantMatrix >= 0.05] = ''

for i, ix in enumerate(plotMatrix.index):
    for j, jx in enumerate(plotMatrix.columns):
        text = g.ax_heatmap.text(
        #text = g.text(
            j + 0.5,
            i + 0.5,
            annot.values[i,j],
            ha="center",
            va="center",
            color="black",
        )
        text.set_fontsize(8)

plt.savefig("results/sigchangefamily.pdf")
#plt.show()
