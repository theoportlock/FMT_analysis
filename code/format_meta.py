#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: Theo Portlock
For diversity analysis
'''
import pandas as pd

# Meta
timemeta = pd.read_csv("data/meta.csv", index_col=0)
meta = timemeta.groupby('id').first()
meta = meta.drop(['MELD','Days after treatment'], axis=1)
meta.to_csv('results/meta.tsv', sep='\t')

var = {'id':'subjectID',
       'Days after treatment':'timepoint',
       'MELD':'MELD'}

timemeta = timemeta[var.keys()]
timemeta = timemeta.rename(columns = var)
timemeta.to_csv('results/timemeta.tsv', sep='\t')
