#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: Theo Portlock
For diversity analysis
'''
import pandas as pd

df = pd.read_csv(f"data/reports/A.Mardinoglu_20_03_sample_info.txt", sep='\t', index_col=0)

ndf = df.Mreads.to_frame()
ndf.index = 'FMT' + ndf.index.str[7:]
ndf.index.name = 'sampleID'

ndf.to_csv(f'results/quality.tsv', sep='\t')
