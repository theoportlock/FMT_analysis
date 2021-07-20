#!/usr/bin/env python
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

# Load data
samples_metadata = pd.read_csv('metadata.csv', index_col=0)
samples_gutmsp = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0).T
samples_oralmsp = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0).T

samples_gutsum = pd.DataFrame(samples_gutmsp.agg(np.count_nonzero, axis=1), columns=['Species richness']).join(samples_metadata)
samples_gutsum['Site'] = 'Gut'
samples_oralsum = pd.DataFrame(samples_oralmsp.agg(np.count_nonzero, axis=1), columns=['Species richness']).join(samples_metadata)
samples_oralsum['Site'] = 'Oral'
joined = pd.concat([samples_gutsum, samples_oralsum]).reset_index()

# Plot
x='Days after treatment'
y='Age'
hue='Type'
col='Site'

g = sns.catplot(data=joined, x=x, y=y, col=col, hue=hue, linewidth=1, kind='box')
#g = sns.catplot(data=joined, x=x, y=y, col=col, hue=hue, linewidth=1, kind='swarm')
#g = sns.catplot(data=joined, x=x, y=y, col=col, hue=hue, linewidth=1, kind='point')

sns.despine(trim=True, left=True)
box_pairs=[((0,'FMT'),(7,'FMT'))]

# Stats
ax, test_results = add_stat_annotation(
    gut_patient_ax,
    data=gut_patient_count_df,
    x=x,
    y=y,
    hue=hue,
    test='Mann-Whitney',
    text_format='full',
    #loc='outside',
    box_pairs=box_pairs,
    verbose=2)

# Print
plt.tight_layout()
plt.show()
#plt.savefig("results/species_richness.pdf")
'''
file_oral_msp_data="../oral_merged_downstream_data/merged.final.mgs.med.vec.10M.csv"
sm = pd.read_csv('metadata.csv')
oral_formatted_names = oral_names['User ID'].str.split(expand=True)
oral_formatted_names.loc[oral_formatted_names[2]=="raw", 2] = 0
oral_formatted_names[2] = oral_formatted_names[2].astype(int)
oral_formatted_names.merge(sm, right_on=['Patient_ID', 'Days after treatment'], left_on=[0,2], how='left')
ofn = oral_formatted_names.merge(sm, right_on=['Patient_ID', 'Days after treatment'], left_on=[0,2], how='left')
oral_msp_data.columns = pd.MultiIndex.from_frame(
        ofn[[0,2,'Type','ID']],
                names=['Patient_ID','Days after treatment', 'Type', 'ID'])
oral_msp_data
oral_msp_data.T
nomd = oral_msp_data.T.reset_index()
fnomd = nomd.iloc[:,3:]
fnomd
ffnomd = fnomd.set_index('ID')
fffnomd = ffnomd.loc[ffnomd.index.dropna()].T
fffnomd.to_csv('oralmsps.csv')
'''
