#!/usr/bin/env python
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from matplotlib_venn import venn2

#samples_metadata = pd.read_csv("../../data/newmergedmetadata.csv").set_index('ID_x')
samples_metadata = pd.read_csv("../../data/newnewmetadata.csv").set_index('Sample ID')
metabolomics = pd.read_csv("../../data/metabolomics.csv",index_col=1).drop('Unnamed: 0', axis=1)

gmsp_samples = pd.read_csv("../../data/gutmsp.csv", index_col=0)
omsp_samples = pd.read_csv("../../data/oralmsp.csv", index_col=0)

gmsp_gtaxonomy = pd.read_csv("../../../downstream_data/taxo.csv", index_col=0)
omsp_otaxonomy = pd.read_csv("../../../oral_downstream_data/oraltaxo.csv", index_col=0)

taxaType='genus'
gtaxonomy_samples = gmsp_samples.join(gmsp_gtaxonomy[taxaType], how='inner').groupby(taxaType).sum()
samples_gtaxonomy = gtaxonomy_samples.T

otaxonomy_samples = omsp_samples.join(omsp_otaxonomy[taxaType], how='inner').groupby(taxaType).sum()
samples_otaxonomy = otaxonomy_samples.T

variables = ['Days after treatment', 'Type', 'Patient_Id'] 
#gjoin = samples_gtaxonomy.join(samples_metadata[variables], how='inner')
#ojoin = samples_otaxonomy.join(samples_metadata[variables], how='inner')
gjoin = samples_gtaxonomy.join(metabolomics, how='inner')
ojoin = samples_otaxonomy.join(metabolomics, how='inner')

gutmetab = gjoin.xs(metabolomics.columns,axis=1)
oralmetab = ojoin.xs(metabolomics.columns,axis=1)
gutnotmetab = gjoin.drop(metabolomics.columns,axis=1)
oralnotmetab = ojoin.drop(metabolomics.columns,axis=1)

c1, p1 = spearmanr(gutmetab, gutnotmetab, axis=0)
c2, p2 = spearmanr(oralmetab, oralnotmetab, axis=0)

'''
ids = metabolomics.merge(samples_metadata, right_on=['time_point', 'Type', 'id'], left_on=['time_point', 'IMP', 'id'])
metabolomics = ids.set_index('ID_x').drop(samples_metadata.drop('ID_x',axis=1).columns, axis=1).drop(['Unnamed: 0', 'Unnamed: 0.1'], axis=1)
metabolomics = ids.copy()

gmerge = gjoin.merge(metabolomics, left_on=variables, right_on=['time_point', 'IMP', 'id'])
'''

correlations1 = pd.DataFrame(
    c1,
    index=gutmetab.columns.append(gutnotmetab.columns),
    columns=gutmetab.columns.append(gutnotmetab.columns))
slicedCorrelations1 = correlations1.iloc[
        len(gutmetab.columns):,
        :len(gutmetab.columns)]

correlations2 = pd.DataFrame(
    c2,
    index=oralmetab.columns.append(oralnotmetab.columns),
    columns=oralmetab.columns.append(oralnotmetab.columns))
slicedCorrelations2 = correlations2.iloc[
        len(oralmetab.columns):,
        :len(oralmetab.columns)]

edges1 = slicedCorrelations1.stack().reset_index()
edges2 = slicedCorrelations2.stack().reset_index()

edges1.columns = ['from','to','value']
edges2.columns = ['from','to','value']

edges1.sort_values('value').tail(20).to_csv('../results/gutmetab.csv',index=False)
edges2.sort_values('value').tail(20).to_csv('../results/oralmetab.csv',index=False)

'''
# for the pvalue number stuff
significantMatrix1 = pd.DataFrame(
    fdrcorrection(slicedCorrelations1.values.flatten())[0].reshape(slicedCorrelations1.shape),
    index = slicedCorrelations1.index,
    columns = slicedCorrelations1.columns)

significantMatrix2 = pd.DataFrame(
    fdrcorrection(slicedCorrelations2.values.flatten())[0].reshape(slicedCorrelations2.shape),
    index = slicedCorrelations2.index,
    columns = slicedCorrelations2.columns)

edges1 = slicedCorrelations1.stack().reset_index()
edges2 = slicedCorrelations2.stack().reset_index()

sigdf = {}
sigdf[0] = {'from':'proteome', 'to':'microbiome', 'spearman':(slicedCorrelations1 > 0.3).sum().sum()}
sigdf[1] = {'from':'proteome', 'to':'metabolome', 'spearman':(slicedCorrelations2 > 0.3).sum().sum()}
sigdf1 = pd.DataFrame(data = sigdf).T

sigdf1.to_csv('total_spearman.csv')

'''
