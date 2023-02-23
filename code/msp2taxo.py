#!/usr/bin/env python
import functions
import pandas as pd
import matplotlib.pyplot as plt

msp = pd.read_csv("../data/gutmsp.csv", index_col=0)
taxo = pd.read_csv("../data/guttaxo.csv", index_col=0)

df = functions.taxofunc(msp, taxo)

df = df.T

meta = pd.read_csv("../data/cleanMeta.csv", index_col=0)

var = ['Days after treatment', 'Type', 'Patient_Id']
#df.join(meta[var]).to_csv('../results/ldadata.tsv', sep='\t')

df1 = df.join(meta[var]).set_index(var)
df1['host_phenotype'] = 'Disease'
disease = df1.xs(0, level=0).set_index('host_phenotype')

df1['host_phenotype'] = 'Healthy'
healthy = df1.xs('DONOR', level=1).set_index('host_phenotype')

df3 = pd.concat([healthy, disease])

df3.T.to_csv("../results/donorldadata.txt", sep="\t")
