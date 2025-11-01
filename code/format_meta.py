#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

# Meta
timemeta = pd.read_csv("data/meta.csv", index_col=0)
meta = timemeta.groupby('id').first()
meta = meta.drop(['MELD','Days after treatment'], axis=1)
meta.index.name = 'subjectID'
meta = meta.rename(columns = {
    'Donor':'donorID',
    'Type':'arm'})
meta['healthy'] = meta.arm == 'DONOR'
meta.to_csv('results/meta.tsv', sep='\t')

var = {'id':'subjectID',
       'Days after treatment':'timepoint',
       'MELD':'MELD'}

timemeta = timemeta[var.keys()]
timemeta = timemeta.rename(columns = var)
timemeta.index.name = 'sampleID'
timemeta.to_csv('results/timemeta.tsv', sep='\t')
