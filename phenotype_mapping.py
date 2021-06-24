#!/usr/bin/env python
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
data = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv",index_col=0)
names = pd.read_table("../downstream_data/P15952_20_03_sample_info_stool.txt")
sample_type = pd.read_csv("../downstream_data/PROFITplaceboTab.csv", index_col=0)
pheno = pd.read_csv("../downstream_data/jgiMat.csv", index_col=0)
variables = pheno.columns
# so that join works
pheno.columns = pd.MultiIndex.from_product([pheno.columns, [""],[""],[""]])

sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))
formatted_names = names['User ID'].str.split(expand=True)
formatted_names.loc[formatted_names[0]=="Donor", 2] = formatted_names.loc[formatted_names[0]=="Donor", 1]
formatted_names.loc[formatted_names[0]=="Donor", 0] = formatted_names.loc[formatted_names[0]=="Donor", 1]
formatted_names.iloc[:,2] = formatted_names.iloc[:,2].str.replace("Stool","0")
formatted_names['ID'] = names['NGI ID'].str.replace('P15952_','FMT')
formatted_names = pd.merge(formatted_names, sample_type, left_on=0, right_on="Patient", how='left')
formatted_names.fillna("DONOR", inplace=True)

data.columns = pd.MultiIndex.from_frame(formatted_names[[0,2,'Type','ID']], names=['Patient','Days after treatment', 'Type', 'ID'])
newpheno=data.join(pheno).dropna()

final = pd.DataFrame(index=data.columns)
final.index = pd.MultiIndex.from_frame(formatted_names[[0,2,'Type','ID']], names=['Patient','Days after treatment', 'Type', 'ID'])
nondonor = final.drop("DONOR", level=2, axis=0)
#donor = final.xs("DONOR", level=2, axis=0)
for i in variables:
    nondonor[i] = newpheno.loc[newpheno[i]==1].sum()
nondonor.index = nondonor.index.set_levels(nondonor.index.levels[1].astype(int), level=1)

groupedfinal = nondonor.groupby(["Days after treatment",'Type']).mean()
groupedfinal = groupedfinal[['Gram+','Gram-','Motile','Sporulating','Nonsporulating','Anaerobe','Facultative','Aerobe','Microaerophilic']]

#groupedfinal.index = groupedfinal.index.astype(int)
groupedfmt = pd.DataFrame(groupedfinal.stack(),columns=["Value"])
groupedfmt.reset_index(inplace=True)
fig = sns.lineplot(
    data=groupedfmt,
    x="Days after treatment",
    y="Value",
    style="Type",
    hue="level_2",
    #palette=sns.color_palette("tab9", len(variables)),
    )
fig.legend(fontsize=8, bbox_to_anchor=(1,1))
plt.tight_layout()
#plt.show()
plt.savefig("results/phenotype_mapping.pdf")

#for i in variables:
#    x,y,hue = "Days after treatment", i, "Type"
#    sns.lineplot(data=nondonor,x=x,y=y,hue=hue)
#    plt.show()
