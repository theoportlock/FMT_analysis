#!/usr/bin/env python
from scipy import stats
from scipy.stats import mannwhitneyu
from scipy.stats import spearmanr
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

samples_metadata = pd.read_csv("../../data/newnewmetadata.csv").set_index('Sample ID')
samples_gutmsp = pd.read_csv("../../data/gutmsp.csv", index_col=0)
samples_oralmsp = pd.read_csv("../../data/oralmsp.csv", index_col=0)
taxaType='genus'
#taxaType='species'
gmsp_taxonomy = pd.read_csv("../../data/guttaxo.csv", index_col=0)
samples_gtaxonomy = samples_gutmsp.join(gmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
gmspdf = samples_gtaxonomy.add_prefix('gut ')
omsp_taxonomy = pd.read_csv("../../data/oraltaxo.csv", index_col=0)
samples_otaxonomy = samples_oralmsp.join(omsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
omspdf = samples_otaxonomy.add_prefix('oral ')
variable=['Days after treatment', 'Type', 'id']
df = omspdf.join(gmspdf, how='inner')
df = df.loc[:, df.nunique() > 10]
df = df.join(samples_metadata[variable], how='inner')
average = df.groupby('Days after treatment').mean()
baseline = df.loc[(df.Type == 'FMT') & (df['Days after treatment'] == 0)].drop("Type", axis =1)
notbaseline  = df.loc[(df.Type == 'FMT') & (df['Days after treatment'] != 0)].drop("Type", axis =1)

# Delta B
v = ['id', 'Days after treatment', 'Type']
taxodata = df.copy()
taxodata = taxodata[taxodata.columns[~taxodata.columns.str.contains('unclassified')]]
tdf = taxodata.reset_index().drop('index', axis=1).set_index(v)
tresult = tdf.div(tdf.reset_index().query('`Days after treatment` == 0').drop(['Days after treatment', 'Type'], axis=1).set_index('id'), level=0).apply(np.log2)
tresult.replace([np.inf, -np.inf], np.nan, inplace=True)
tresult.fillna(0, inplace=True)
tmeanresult = tresult.reset_index().drop('id', axis=1).groupby(['Days after treatment', 'Type']).mean()
tfiltmeanresult = tmeanresult.loc[:, (np.abs(stats.zscore(tmeanresult)) > 2.6).any()]
tfiltmeanresult.drop(0, level=0, inplace=True)
sns.heatmap(tfiltmeanresult.T, center=0 ,cmap='vlag')
tdf.reset_index(inplace=True)
tdf.drop('id', axis=1, inplace=True)
variables=['Days after treatment', 'Type'] 
baseline = tdf.loc[(tdf['Type'] == 'FMT') & (tdf['Days after treatment'] == 0)].drop(variables, axis =1)
#baseline = tdf.loc[tdf['Days after treatment'] == 0].drop('Days after treatment', axis=1).set_index('Type')
notbaseline = tdf.loc[(tdf['Type'] == 'FMT') & (tdf['Days after treatment'] != 0)].drop('Type', axis=1)
#notbaseline = tdf.loc[tdf['Days after treatment'] != 0].set_index(['Type', 'Days after treatment'])
def mandf(df):
    base = baseline
    day = df
    result = pd.Series()
    for i in base.columns:
        try:
            result[i] = mannwhitneyu(base[i],day[i])[1]
        except:
            result[i] = np.nan
    return result

#notbaseline.groupby(level=[0, 1]).apply(lambda t: stats.mannwhitneyu(baseline.loc[t], t)[1][1])
pvals = pd.DataFrame(notbaseline.groupby('Days after treatment').apply(mandf))
qvals = pd.DataFrame(fdrcorrection(pvals.values.flatten())[0].reshape(pvals.shape), index = pvals.index, columns = pvals.columns)
bsignificantMatrix = pvals.loc[:, (pvals < 0.05).any(axis=0)]
#significantMatrix = qvals.loc[:, (qvals < 0.05).any(axis=0)]
plotdf = tmeanresult[bsignificantMatrix.columns]
plotdf.drop(0, level=0, inplace=True)
sns.heatmap(plotdf.T, cmap='vlag')

# Delta P
v = ['id', 'Days after treatment', 'Type']
taxodata = df.copy()
taxodata = taxodata[taxodata.columns[~taxodata.columns.str.contains('unclassified')]]
tdf = taxodata.reset_index().drop('index', axis=1).set_index(v)
tresult = tdf.div(tdf.reset_index().query('`Days after treatment` == 0').drop(['Days after treatment', 'Type'], axis=1).set_index('id'), level=0).apply(np.log2)
tresult.replace([np.inf, -np.inf], np.nan, inplace=True)
tresult.fillna(0, inplace=True)
tmeanresult = tresult.reset_index().drop('id', axis=1).groupby(['Days after treatment', 'Type']).mean()
tfiltmeanresult = tmeanresult.loc[:, (np.abs(stats.zscore(tmeanresult)) > 2.6).any()]
tfiltmeanresult.drop(0, level=0, inplace=True)
sns.heatmap(tfiltmeanresult.T, center=0 ,cmap='vlag')
tdf.reset_index(inplace=True)
tdf.drop('id', axis=1, inplace=True)
variables=['Days after treatment', 'Type'] 
baseline = tdf.loc[tdf['Type'] == 'FMT'].drop('Type', axis=1)
#baseline = tdf.loc[tdf['Days after treatment'] == 0].drop('Days after treatment', axis=1).set_index('Type')
notbaseline = tdf.loc[tdf['Type'] != 'FMT'].drop('Type', axis=1)
#notbaseline = tdf.loc[tdf['Days after treatment'] != 0].set_index(['Type', 'Days after treatment'])
def nmandf(df):
    var = df['Days after treatment'].unique()[0]
    base = baseline.drop('Days after treatment', axis=1)
    day = df
    result = pd.Series()
    for i in base.columns:
        try:
            result[i] = mannwhitneyu(base[i],day[i])[1]
        except:
            result[i] = np.nan
    return result

#notbaseline.groupby(level=[0, 1]).apply(lambda t: stats.mannwhitneyu(baseline.loc[t], t)[1][1])
notbaseline.groupby(level=[0, 1]).apply(lambda t: stats.mannwhitneyu(baseline.loc[t], t)[1][1])
pvals = pd.DataFrame(notbaseline.groupby('Days after treatment').apply(nmandf))
qvals = pd.DataFrame(fdrcorrection(pvals.values.flatten())[0].reshape(pvals.shape), index = pvals.index, columns = pvals.columns)
psignificantMatrix = pvals.loc[:, (pvals < 0.05).any(axis=0)]
#significantMatrix = qvals.loc[:, (qvals < 0.05).any(axis=0)]
plotdf = tmeanresult[significantMatrix.columns]
plotdf.drop(0, level=0, inplace=True)

presult = tdf.reset_index().drop('id', axis=1).groupby(['Days after treatment', 'Type']).mean()
ppresult = presult.xs('FMT', level=1).div(presult.xs('PLACEBO', level=1)).apply(np.log2)
ppresult.replace([np.inf, -np.inf], np.nan, inplace=True)
ppresult.fillna(0, inplace=True)
#pmeanresult = presult.reset_index().drop('id', axis=1).groupby(['Days after treatment', 'Type']).mean()
plotdf = ppresult[psignificantMatrix.columns]
sns.heatmap(plotdf.T, cmap='vlag', yticklabels=True)

# Delta D
v = ['id', 'Days after treatment', 'Type']
taxodata = df.copy()
taxodata = taxodata[taxodata.columns[~taxodata.columns.str.contains('unclassified')]]
tdf = taxodata.reset_index().drop('index', axis=1).set_index(v)
tresult = tdf.div(tdf.reset_index().query('`Days after treatment` == 0').drop(['Days after treatment', 'Type'], axis=1).set_index('id'), level=0).apply(np.log2)
tresult.replace([np.inf, -np.inf], np.nan, inplace=True)
tresult.fillna(0, inplace=True)
tmeanresult = tresult.reset_index().drop('id', axis=1).groupby(['Days after treatment', 'Type']).mean()
tfiltmeanresult = tmeanresult.loc[:, (np.abs(stats.zscore(tmeanresult)) > 2.6).any()]
tfiltmeanresult.drop(0, level=0, inplace=True)
sns.heatmap(tfiltmeanresult.T, center=0 ,cmap='vlag')
tdf.reset_index(inplace=True)
tdf.drop('id', axis=1, inplace=True)
variables=['Days after treatment', 'Type'] 
baseline = gmspdf.join(samples_metadata[variables], how='inner')
baseline = baseline.loc[baseline['Type'] == 'DONOR'].drop(variables, axis=1)
notbaseline = gmspdf.join(samples_metadata[variables], how='inner')
notbaseline = notbaseline.loc[notbaseline['Type'] != 'DONOR']
#notbaseline = tdf.loc[tdf['Type'] != 'FMT'].drop('Type', axis=1)
#notbaseline = tdf.loc[tdf['Days after treatment'] != 0].set_index(['Type', 'Days after treatment'])
def dmandf(df):
    var = df['Days after treatment'].unique()[0]
    base = baseline
    day = df
    result = pd.Series()
    for i in base.columns:
        try:
            result[i] = mannwhitneyu(base[i],day[i])[1]
        except:
            result[i] = np.nan
    return result

notbaseline.groupby(level=[0, 1]).apply(lambda t: stats.mannwhitneyu(baseline.loc[t], t)[1][1])
pvals = pd.DataFrame(notbaseline.groupby(variables).apply(dmandf))
#qvals = pd.DataFrame(fdrcorrection(pvals.values.flatten())[0].reshape(pvals.shape), index = pvals.index, columns = pvals.columns)
dsignificantMatrix = pvals.loc[:, (pvals < 0.05).any(axis=0)]
#significantMatrix = qvals.loc[:, (qvals < 0.05).any(axis=0)]

dresult = tdf.groupby(['Days after treatment', 'Type']).mean()
ddresult = dresult.div(baseline.mean()).apply(np.log2)
ddresult.replace([np.inf, -np.inf], np.nan, inplace=True)
ddresult.fillna(0, inplace=True)
#pmeanresult = presult.reset_index().drop('id', axis=1).groupby(['Days after treatment', 'Type']).mean()
plotdf = ddresult[dsignificantMatrix.columns]
sns.heatmap(plotdf.T, cmap='vlag', yticklabels=True)
