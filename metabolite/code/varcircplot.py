#!/usr/bin/env python
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy import stats

#metab = pd.read_csv("../data/serum_bileacid.csv", index_col=2)
metab = pd.read_csv("../data/faecal_bileacids.csv", index_col=0)
meta = pd.read_csv("../data/carbalivemeta.tsv", sep='\t').dropna(subset=['Gnumber']).set_index('Gnumber')
data = pd.read_csv("../data/msp.csv", index_col=0).T
taxo = pd.read_csv("../data/taxo.csv", index_col=0)
data.index = data.index.str.replace(r'.*_G','G').str.replace(r'_*','')
data.index = data.index.str.replace(r'S.*', '')

metab.index = metab.index.str[:-2]
#metab.drop(['Unnamed: 0', 'Unnamed: 1', 'Campione'], axis=1, inplace=True)
metab.drop(['ID'], axis=1, inplace=True)
metab = metab.loc[:, metab.sum() != 0]
taxa_type='species'
variables = ['SUBJID','VISIT', 'ARM']
df = metab.join(meta[variables] , how='inner').set_index(variables)

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
taxodata = data.T.join(taxo[taxa_type], how='right').groupby(taxa_type).sum().T.fillna(0)
taxodata = taxodata[taxodata.columns[~taxodata.columns.str.contains('unclassified')]]
taxodata = taxodata.loc[:, taxodata.nunique() > 10]
tdf = taxodata.join(meta[variables], how='inner').set_index(variables)
tresult = tdf.div(tdf.reset_index().query('VISIT == "BASELINE"').drop(['VISIT', 'ARM'], axis=1).set_index('SUBJID'), level=0).apply(np.log2)
tresult.replace([np.inf, -np.inf], np.nan, inplace=True)
tresult.fillna(0, inplace=True)
tmeanresult = tresult.reset_index().drop('SUBJID', axis=1).groupby(['VISIT', 'ARM']).mean()
tfiltmeanresult = tmeanresult.loc[:, (np.abs(stats.zscore(tmeanresult)) > 2.15).any()]
#tfiltmeanresult = tmeanresult.copy()
tfiltmeanresult.drop('BASELINE', level=0, inplace=True)
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
tdf.drop('SUBJID', axis=1, inplace=True)
variables=['VISIT', 'ARM'] 
#tdf = tdf.reset_index()
baseline = tdf.loc[(tdf.ARM == 'ACTIVE') & (tdf['VISIT'] == 'BASELINE')].drop(variables, axis =1)
notbaseline = tdf.loc[(tdf.ARM == 'ACTIVE') & (tdf['VISIT'] != 'BASELINE')].drop('ARM', axis=1)
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
stats = pd.DataFrame(notbaseline.groupby('VISIT').apply(mandf))
#classstats = stats[stats.columns[~stats.columns.str.contains('unclassified')]].drop('index', axis=1)
#classstats = stats.loc[:, stats.dtypes == np.float64].drop('index', axis=1)
                
significantMatrix = stats.loc[:, (stats < 0.15).any(axis=0)]
plotdf = tmeanresult[significantMatrix.columns]
plotdf.drop('BASELINE', level=0, inplace=True)
sns.heatmap(plotdf.T, cmap='vlag')
