#!/usr/bin/env python
import gc
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import dask.dataframe as dd

from dask.distributed import Client
# set up cluster and workers
client = Client(n_workers=4, threads_per_worker=1, memory_limit='64GB') 
client

# Data location
file_gut_msp_data="../project/Downstream/norm.csv"
file_gene_info='../downstream_data/IGC2.1990MSPs.tsv'
file_kegg='../downstream_data/hs_10_4_igc2_annot_20180614_vs_kegg82.best.csv'
file_patric='../downstream_data/patricmap.csv'
file_antismash='../downstream_data/igc2.antismash.simple.txt'
file_pfam='../downstream_data/igc2PfamDesc.csv'
file_patric='../downstream_data/patricmap.csv'
file_card='../downstream_data/hs_10_4_igc2.CARD.tsv'

# Load data
gene_info = dd.read_csv(file_gene_info, sep='\t')
gut_msp_data = dd.read_csv(file_gut_msp_data)

# Map to catalogues
## KEGG
kegg = dd.read_csv(file_kegg)
kegg.columns = ['gene_id','gene_name','ko', 'bitscore']
replacedbykegg = gut_msp_data.merge(kegg, how='inner', left_on='Unnamed: 0', right_on='gene_id')
fmtcols = [col for col in replacedbykegg.columns.values if 'FMT' in col]
fmtcols.append('ko')
sumofsm = replacedbykegg[fmtcols].groupby('ko').sum()
del replacedbykegg
del kegg
gc.collect()
sumofsm.to_csv("keggNorm.csv",  single_file=True)
del sumofsm


## antismash
antismash = pd.read_table(file_antismash)
antismash.columns = ['gene_name','sm','id']
replacedbygenename = gut_msp_data.merge(gene_info, how='inner', left_on='Unnamed: 0', right_on='gene_id')
replacedbyantismash = replacedbygenename.merge(antismash, how='inner', left_on='gene_name', right_on='gene_name')
fmtcols = [col for col in replacedbyantismash.columns.values if 'FMT' in col]
fmtcols.append('sm')
sumofsm = replacedbyantismash[fmtcols].groupby('sm').sum()
del replacedbyantismash, antismash
gc.collect()
sumofsm.to_csv("antismashNorm.csv",  single_file=True)
del sumofsm


## pfam 
pfam = dd.read_csv(file_pfam, assume_missing=True)
replacedbygenename = gut_msp_data.merge(gene_info, how='inner', left_on='Unnamed: 0', right_on='gene_id')
replacedbypfam = replacedbygenename.merge(pfam, how='inner', left_on='gene_name', right_on='gene_name')
fmtcols = [col for col in replacedbypfam.columns.values if 'FMT' in col]
fmtcols.append('pfam_name')
sumofsm = replacedbypfam[fmtcols].groupby('pfam_name').sum()
del replacedbypfam, pfam
gc.collect()
sumofsm.to_csv("pfamNorm.csv",  single_file=True)
del sumofsm


## patric 
patric = dd.read_csv(file_patric)
replacedbygenename = gut_msp_data.merge(gene_info, how='inner', left_on='Unnamed: 0', right_on='gene_id')
replacedbypatric = replacedbygenename.merge(patric, how='inner', left_on='gene_name', right_on='igc2_id')
del replacedbygenename
fmtcols = [col for col in replacedbypatric.columns.values if 'FMT' in col]
fmtcols.append('vf_id')
sumofsm = replacedbypatric[fmtcols].groupby('vf_id').sum()
del replacedbypatric, patric
gc.collect()
sumofsm.to_csv("patricNorm.csv",  single_file=True)
del sumofsm


## CARD
card = pd.read_csv(file_card, sep='\t')
card["gene_name"] = card.ORF_ID.str.split(expand=True)[0]
card["gene_name"] = card["gene_name"].str[0:-2]
replacedbygenename = gut_msp_data.merge(gene_info, how='inner', left_on='Unnamed: 0', right_on='gene_id')
replacedbycard = replacedbygenename.merge(card, how='inner', left_on='gene_name', right_on='gene_name')
del replacedbygenename
fmtcols = [col for col in replacedbycard.columns.values if 'FMT' in col]
fmtcols.append('Best_Hit_ARO')
sumofsm = replacedbycard[fmtcols].groupby('Best_Hit_ARO').sum()
del replacedbycard, card
gc.collect()
sumofsm.to_csv("cardNorm.csv",  single_file=True)
del sumofsm

client.close()
