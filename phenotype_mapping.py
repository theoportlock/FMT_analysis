#!/usr/bin/env python
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

sns.set(rc={'figure.figsize':(11.7,8.27)})

taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
data = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv",index_col=0)
pheno = pd.read_csv("../downstream_data/jgiMat.csv", index_col=0)

variables = pheno.columns
samples_metadata = pd.read_csv('metadata.csv', index_col=0)

newpheno=data.join(pheno, how='inner').dropna()

nondonor = pd.DataFrame(index=data.columns)

for i in variables: nondonor[i] = newpheno.loc[newpheno[i]==1].sum()

final = nondonor.join(samples_metadata[["Days after treatment", "Type"]], how='inner')

plotfinal = final.loc[final.Type == 'FMT'].drop(['Type'],axis=1).groupby("Days after treatment").mean().T

sns.heatmap(plotfinal, cmap='binary')

plt.tight_layout()
plt.show()
plt.savefig("results/new_phenotype_mapping.pdf")
