#!/usr/bin/env python
from scipy import stats
from scipy.stats import spearmanr
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# df preparation
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

# BA FOLD CHANGE
result = df.div(df.reset_index().query('VISIT == "BASELINE"').drop(['VISIT', 'ARM'], axis=1).set_index('SUBJID'), level=0).apply(np.log2)
result.replace([np.inf, -np.inf], np.nan, inplace=True)
result.fillna(0, inplace=True)
meanresult = result.reset_index().drop('SUBJID', axis=1).groupby(['VISIT', 'ARM']).mean()
#filtmeanresult = meanresult.loc[:, (np.abs(stats.zscore(meanresult)) > 2).any()]
filtmeanresult = meanresult.copy()
filtmeanresult.drop('BASELINE', level=0, inplace=True)
#plt.rcParams["figure.figsize"] = (4,4)
#sns.clustermap(filtmeanresult.T, center=0 ,cmap='vlag')
sns.heatmap(filtmeanresult.T, center=0 ,cmap='vlag')

# SPEC FOLD CHANGE - deffo correct
#taxodata = data.T.join(taxo[taxa_type], how='right').groupby(taxa_type).sum().T.fillna(0)
v = ['id', 'Days after treatment', 'Type']
taxodata = df.copy()
taxodata = taxodata[taxodata.columns[~taxodata.columns.str.contains('unclassified')]]
#tdf = taxodata.join(meta[variables], how='inner').set_index(variables)
tdf = taxodata.reset_index().drop('index', axis=1).set_index(v)

#tresult = tdf.div(tdf.reset_index().query('`Days after treatment` == 0').drop(["Days after treatment", 'Type'], axis=1).set_index('SUBJID'), level=0).apply(np.log2)
#tresult = tdf.div(tdf.reset_index().query('`Days after treatment` == 0').drop(v, axis=1), level=0).apply(np.log2)
tresult = tdf.div(tdf.reset_index().query('`Days after treatment` == 0').drop(['Days after treatment', 'Type'], axis=1).set_index('id'), level=0).apply(np.log2)
#tresult = tdf.div(tdf.query('`Days after treatment` == 0', level=1).apply(np.log2)
tresult.replace([np.inf, -np.inf], np.nan, inplace=True)
tresult.fillna(0, inplace=True)
tmeanresult = tresult.reset_index().drop('id', axis=1).groupby(['Days after treatment', 'Type']).mean()
tfiltmeanresult = tmeanresult.loc[:, (np.abs(stats.zscore(tmeanresult)) > 2.6).any()]
#tfiltmeanresult = tmeanresult.copy()
tfiltmeanresult.drop(0, level=0, inplace=True)
sns.heatmap(tfiltmeanresult.T, center=0 ,cmap='vlag')

# Corrs
sdata = data.T.join(taxo['species']).groupby('species').sum().T
sdata = sdata.loc[:, sdata.nunique() > 10]
join = sdata.join(metab, how='inner')
join = join.loc[:, ~join.columns.str.contains('unclassified')]
fmetab = join.xs(metab.columns, axis=1)
fnotmetab = join.drop(metab.columns, axis=1)
c, p = spearmanr(join, axis=0)
corrs = join.corr(method='spearman')
#corrs = pd.DataFrame(
    #c,
    #index=fmetab.columns.append(fnotmetab.columns),
    #columns=fmetab.columns.append(fnotmetab.columns))
slicedcorrs = corrs.iloc[
    len(fnotmetab.columns):,
    :len(fnotmetab.columns)]
edges = slicedcorrs.stack().reset_index()
edges.columns = ['from','to','value']

slicededges = edges.loc[np.abs(edges.value) > 0.4]
slicededges.to_csv('snewedges.csv', index=False)
#tfiltmeanresult.T['WEEK03', 'ACTIVE'].to_csv('snodes.csv')
#tf = tfiltmeanresult.T.loc[slicededges['from'].unique()]['WEEK04', 'ACTIVE']
tf = tfiltmeanresult.T.loc[slicededges['from'].unique()]
#filtmeanresult.T['WEEK04', 'ACTIVE'].to_csv('mnodes.csv')
#fi = filtmeanresult.T.loc[slicededges['to'].unique()]['WEEK04', 'ACTIVE']
fi = filtmeanresult.T.loc[slicededges['to'].unique()]
tf.append(fi).to_csv('nnodes.csv')

# Slicedheatmap
#sns.clustermap(slicedcorrs[tfiltmeanresult.columns], cmap='vlag')
#sns.clustermap(slicedcorrs[significantMatrix.columns], cmap='vlag')
sslicedcorrs = slicedcorrs[significantMatrix.columns]
sedges = sslicedcorrs.stack().reset_index()
sedges.columns = ['from','to','value']
#sedges.to_csv('ffnewedges.csv', index=False)
sedges[np.abs(sedges.value) > 0.2].to_csv('ffnewedges.csv', index=False)

tmp = corrs.iloc[:len(fnotmetab.columns),:len(fnotmetab.columns)]
tedges = tmp.stack().reset_index()
tedges.columns = ['from','to','value']
tedges[np.abs(tedges.value) > 0.5].to_csv('tedges.csv', index=False)

# Sigtaxachange
from scipy.stats import mannwhitneyu
tdf.reset_index(inplace=True)
tdf.drop('id', axis=1, inplace=True)
variables=['Days after treatment', 'Type'] 
#tdf = tdf.reset_index()
baseline = tdf.loc[(tdf['Type'] == 'FMT') & (tdf['Days after treatment'] == 0)].drop(variables, axis =1)
notbaseline = tdf.loc[(tdf['Type'] == 'FMT') & (tdf['Days after treatment'] != 0)].drop('Type', axis=1)
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
pvals = pd.DataFrame(notbaseline.groupby('Days after treatment').apply(mandf))
#classpvals = pvals[pvals.columns[~pvals.columns.str.contains('unclassified')]].drop('index', axis=1)
#classpvals = pvals.loc[:, pvals.dtypes == np.float64].drop('index', axis=1)
                
qvals = pd.DataFrame(
    fdrcorrection(pvals.values.flatten())[0].reshape(pvals.shape),
    index = pvals.index,
    columns = pvals.columns)

significantMatrix = qvals.loc[:, (qvals < 0.05).any(axis=0)]

plotdf = tmeanresult[significantMatrix.columns]
plotdf.drop(0, level=0, inplace=True)
sns.heatmap(plotdf.T, cmap='vlag')
