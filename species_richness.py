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
joined['Type'] = joined.Type.fillna('DONOR')
joined.loc[joined.Type == 'DONOR', 'Site'] = 'DONOR'

# Plot
x='Days after treatment'
y='Species richness'
hue='Type'
col='Site'

#g = sns.catplot(data=joined, x=x, y=y, col=col, hue=hue, linewidth=1, kind='box')
fig, (gut_donor_ax, gut_patient_ax, oral_patient_ax) = plt.subplots(1, 3, sharey=True, gridspec_kw={'width_ratios': [1, 7, 7]})
sns.set_theme(font_scale=0.8)
sns.boxplot(ax = gut_donor_ax, data=joined.loc[joined.Site == 'DONOR'], x="Type", y=y, linewidth=1).set_title("Donor stool")
sns.stripplot(ax = gut_donor_ax, data=joined.loc[joined.Site == 'DONOR'], x="Type", y=y, linewidth=1)
sns.boxplot(ax = gut_patient_ax, data=joined.loc[joined.Site == 'Gut'], x=x, y=y, hue=hue, linewidth=1).set_title("Patient stool")
sns.stripplot(ax = gut_patient_ax, data=joined.loc[joined.Site == 'Gut'], x=x, y=y, hue=hue, linewidth=1)
sns.boxplot(ax = oral_patient_ax, data=joined.loc[joined.Site == 'Oral'], x=x, y=y, hue=hue, linewidth=1).set_title("Patient Saliva")
sns.stripplot(ax = oral_patient_ax, data=joined.loc[joined.Site == 'Oral'], x=x, y=y, hue=hue, linewidth=1)

# Format Plot
gut_patient_ax.legend([],[], frameon=False)
gut_donor_ax.legend([],[], frameon=False)
gut_patient_ax.set_ylabel('')
oral_patient_ax.set_ylabel('')
gut_donor_ax.set_xlabel('')
sns.despine(trim=True, left=True)
box_pairs=[((0,'FMT'),(7,'FMT'))]

# Stats
ax, test_results = add_stat_annotation(
    gut_patient_ax,
    data=joined.loc[joined.Site == 'Gut'],
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
