# coding: utf-8
#!/usr/bin/env python
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import math
from bioinfokit import analys, visuz
from scipy.stats import wilcoxon as wc
 
pd.set_option('use_inf_as_na', True)

samples_metadata = pd.read_csv('metadata.csv', index_col=0).dropna()
samples_patric = pd.read_csv("patricNorm.csv", index_col=0).T
df = samples_patric.join(samples_metadata[["Days after treatment", "Type"]], how='inner')

#df = pickle.load(open("../downstream/Func.pickle",'rb'))
#df2 = df.reset_index().copy()

#df3 = df2.drop("ID",axis="columns")
#df3.set_index(['Patient','Days after treatment', 'Type'], inplace=True)

#df3.drop(["P0017","P0023","P0029"], inplace=True)
    
#df3 = df3[df3.columns[df3.sum()>0.001]].copy()
day=7
Fd7='`Days after treatment` == @day and Type == "FMT"'
Fd0='`Days after treatment` == 0 and Type == "FMT"'
df.query(Fd0)
df.query(Fd0).apply(np.log)

logdif = df.query(Fd0).apply(np.log) - df.query(Fd0).apply(np.log)
logdif = df3.xs([0,'FMT'],axis='index',level=[1,2]).apply(np.log) - df3.xs([day,'FMT'],axis='index',level=[1,2]).apply(np.log)

(df.xs([0,'FMT']).apply(np.log).mean() - df.xs([7,'FMT']).apply(np.log).mean()).dropna()

comb = pd.DataFrame()
comb['logFC'] = logdif.mean()

minlogtenpval = pd.DataFrame()
for i in df3.columns:
    minlogtenpval[i] = -np.log(wc(df3.xs([0,'FMT'],axis='index',level=[1,2])[i], df3.xs([day,'FMT'],axis='index',level=[1,2])[i]))
    
    
comb['-log10(Pvalue)'] = minlogtenpval.T[1].copy()

comb.replace([np.inf, -np.inf], np.nan, inplace=True)
comb.dropna(axis='index',inplace=True)

pval = pd.DataFrame()
for i in df3.columns:
    pval[i] = wc(df3.xs([0,'FMT'],axis='index',level=[1,2])[i], df3.xs([day,'FMT'],axis='index',level=[1,2])[i])
pval

comb['P-value'] = pval.T[1].copy()
visuz.gene_exp.volcano(df=comb.reset_index(), lfc='logFC', pv='P-value', geneid='index', sign_line=True, show=True, plotlegend=True, legendpos='upper right', genenames='deg')
