#!/usr/bin/env python
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import math
from bioinfokit import analys, visuz
from scipy.stats import wilcoxon as wc

df = pickle.load(open("data/Func.pickle",'rb'))
df2 = df.reset_index().copy()

df3 = df2.drop("ID",axis="columns")
df3.set_index(['Patient','Days after treatment', 'Type'], inplace=True)

df3.drop(["P0017","P0023","P0029"], inplace=True)
    
df3 = df3[df3.columns[df3.sum()>0.001]].copy()

comb = pd.DataFrame()

days=[7,30,90]

for day in days:
    logdif = df3.xs([0,'FMT'],axis='index',level=[1,2]).apply(np.log) - df3.xs([day,'FMT'],axis='index',level=[1,2]).apply(np.log)

    #logdif = df3.xs([0,'FMT'],axis='index',level=[1,2]) - df3.xs([day,'FMT'],axis='index',level=[1,2])

    comb['logFC'] = logdif.mean()
    comb.replace([np.inf, -np.inf], np.nan, inplace=True)
    #comb.dropna(axis='index',inplace=True)

    pval = pd.DataFrame()
    for i in df3.columns:
        pval[i] = wc(df3.xs([0,'FMT'],axis='index',level=[1,2])[i], df3.xs([day,'FMT'],axis='index',level=[1,2])[i])

    comb['P-value'] = pval.T[1].copy()
