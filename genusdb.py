#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Plotting an annotated seaborn correlation heatmap between MSP data and metadata
Theo Portlock
%autoindent
'''

import pandas as pd
msp_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
msp_samples = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)

# Join taxo information
taxaType='genus'
taxonomy_samples = msp_samples.join(msp_taxonomy[taxaType], how='inner').groupby(taxaType).sum()
taxonomy_samples.to_csv('genusdb.tsv')
