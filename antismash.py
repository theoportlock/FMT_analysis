#!/usr/bin/env python
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import dask.dataframe as dd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

# Data location
file_gut_taxonomy="../downstream_data/taxo.csv"
file_gut_msp_data="../downstream_data/gctDown10m.csv"
file_gut_names="../downstream_data/P15952_20_03_sample_info_stool.txt"
file_sample_type="../downstream_data/PROFITplaceboTab.csv"
file_antismash='../downstream_data/igc2.antismash.simple.txt'
file_hs='../downstream_data/hs_v2_sp.csv'

# Load data
sample_type = pd.read_csv(file_sample_type, index_col=0)

gut_taxonomy = pd.read_csv(file_gut_taxonomy, index_col=0)
gut_msp_data = dd.read_csv(file_gut_msp_data)
gut_names = pd.read_table(file_gut_names, index_col=0)
antismash = pd.read_table(file_antismash)
antismash.columns = ['gene_name','sm','id']
antismash.set_index('gene_name', inplace=True)
hs = pd.read_csv(file_hs, index_col = 0)

sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))

gut_names["User ID"] = gut_names["User ID"].str.replace("Baseline", "Day 0")
gut_formatted_names = gut_names['User ID'].str.split(expand=True)
gut_formatted_names.loc[gut_formatted_names[0]=="Donor", 2] = gut_formatted_names.loc[gut_formatted_names[0]=="Donor", 1]
gut_formatted_names.iloc[:,2] = gut_formatted_names.iloc[:,2].str.replace("Stool","0")
gut_formatted_names = pd.merge(gut_formatted_names, sample_type, left_on=0, right_on="Patient", how='left')
gut_msp_data = gut_msp_data.set_index("Unnamed: 0")
gut_msp_data.columns = pd.MultiIndex.from_frame(
        gut_formatted_names[[0,2,'Type']],
        names=['Patient','Days after treatment', 'Type'])

dtaxa_not_sum = gut_msp_data.join(gut_taxonomy['gene_name']).set_index('gene_name')
replacedbygenename = gut_msp_data.join(hs['gene_name'])
