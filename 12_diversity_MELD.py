#!/usr/bin/env python
import os
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation
from skbio.diversity.alpha import shannon

file_taxonomy="../downstream_data/taxo.csv"
file_meteor_data="../downstream_data/merged.final.mgs.med.vec.10M.csv"
file_names="../downstream_data/P15952_20_03_sample_info_stool.txt"
file_sample_type="../downstream_data/PROFITplaceboTab.csv"
file_day0="../downstream_data/PROFITbaseline.csv"
file_day7="../downstream_data/PROFITday7.csv"
file_day30="../downstream_data/PROFITday30.csv"
file_day90="../downstream_data/PROFITday90.csv"
file_siafile="../downstream_data/normalised.faecal.cytokines.markers.xlsx"

taxonomy = pd.read_csv(file_taxonomy, index_col=0)
data = pd.read_csv(file_meteor_data, index_col=0)
names = pd.read_table(file_names, index_col=0)
sample_type = pd.read_csv(file_sample_type, index_col=0)
day0=pd.read_csv(file_day0)
day7=pd.read_csv(file_day7)
day30=pd.read_csv(file_day30)
day90=pd.read_csv(file_day90)
siafile = pd.read_excel(file_siafile)

sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))
formatted_names = names['User ID'].str.split(expand=True)
formatted_names.loc[formatted_names[0]=="Donor", 2] = formatted_names.loc[formatted_names[0]=="Donor", 1]
formatted_names.iloc[:,2] = formatted_names.iloc[:,2].str.replace("Stool","0")
formatted_names = pd.merge(formatted_names, sample_type, left_on=0, right_on="Patient", how='left')
data.columns = pd.MultiIndex.from_frame(
        formatted_names[[0,2,'Type']],
        names=['Patient','Days after treatment', 'Type'])
data = data.drop("Donor", axis=1, level=0)
count_df = pd.DataFrame(data.agg(shannon, axis=0)).reset_index()
count_df = count_df.rename(columns={0:"Shannon diversity"})
count_df["Days after treatment"] = count_df["Days after treatment"].astype(int)

day0["Days after treatment"]=0
day7["Days after treatment"]=7
day30["Days after treatment"]=30
day90["Days after treatment"]=90
vals = pd.concat([day0,day7,day30,day90])
vals["Patient_Id"] = vals.Patient.str.replace("P01","P").str.replace("\ .*", "", regex=True)
siafile["id"] = siafile.id.str.replace("P01","P").str.replace("\ .*", "", regex=True)
siafile["time_point"] = siafile.time_point.str.replace("D","")
siafile['time_point'] = siafile.time_point.astype('int64')
merged_metadata = pd.merge(siafile, vals,  how='inner', left_on=['id','time_point'], right_on = ['Patient_Id','Days after treatment'])

merged_metadata_diversity = pd.merge(merged_metadata, count_df, how='inner', left_on=['id','time_point'], right_on=['Patient', 'Days after treatment'])

sns.regplot(data=merged_metadata_diversity,x='MELD_y', y='Shannon diversity')

'''
x='Days after treatment'
y='Shannon diversity'
hue='Type'

ax = sns.boxplot(data=count_df,x=x,y=y, hue=hue, linewidth=1)

box_pairs=[((0,'FMT'),(7,'FMT'))]
#box_pairs = list(itertools.combinations([(a,b) for a in count_df['Days after treatment'].unique() for b in count_df["Type"].unique()],2))

ax, test_results = add_stat_annotation(
    ax,
    data=count_df,
    x=x,
    y=y,
    hue=hue,
    test='Mann-Whitney',
    #test='Wilcoxon',
    text_format='full',
    #loc='outside',
    box_pairs=box_pairs,
    verbose=2)

sns.stripplot(data=count_df,x=x,y=y, hue=hue, linewidth=1,size=4)
sns.despine(trim=True, left=True)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
ax.set_title("Change in gut microbiome species diversity after FMT treatment")
plt.tight_layout()

plt.savefig("results/1_shannon_diversity.pdf")
'''
plt.show()
#plt.savefig("results/" + os.path.basename(__file__) + ".pdf")
