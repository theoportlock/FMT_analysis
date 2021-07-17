#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
A script for plotting an annotated seaborn heatmap
Theo Portlock
'''
%autoindent

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import plot_correlation.main as plotcor

metadata = pd.read_csv('metadata.csv', index_col=0)
taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
msp_data = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
gut_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)

# Join taxo information
taxa_type='species'
taxa_sampleid = msp_data.join(gut_taxonomy[taxa_type])
msp_data = taxa_sampleid.groupby(taxa_type).sum().T
merged_msp_metadata = msp_data.join(metadata)
numeric_data = merged_msp_metadata.dropna()._get_numeric_data()

plotcor(numeric_data)

plt.show()
plt.tight_layout()
plt.savefig("clustermap.png", dpi=200)
