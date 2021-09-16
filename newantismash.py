#!/usr/bin/env python
import gc
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import dask.dataframe as dd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

samples_metadata = pd.read_csv('newmergedmetadata.csv').drop_duplicates(subset='ID_x', keep='first').set_index('ID_x').sort_index().dropna(thresh=20, axis=1)
samples_sm = pd.read_csv("antismashNorm.csv", index_col=0).T

samples_metadata.loc[samples_metadata.Type == 'DONOR', 'Days after treatment'] = 0

samples_smMetadata = samples_sm.join(samples_metadata[['Days after treatment', 'Type']], how='inner')
samples_smMetadata = samples_smMetadata.groupby(['Days after treatment', 'Type']).mean()
samples_smMetadata = samples_smMetadata/samples_smMetadata.max().max()

# plot
sns.set(rc={'figure.figsize':(11.7,8.27)})
fig, (gut_donor_ax, gut_patient_fmt_ax, gut_patient_placebo_ax) = plt.subplots(1, 3, sharey=True, gridspec_kw={'width_ratios': [1, 4, 4]})

# donor
datt3 = samples_smMetadata.query('Type == "DONOR"').droplevel('Type').T
sns.heatmap(datt3,cmap='binary',ax=gut_donor_ax, cbar=False)
gut_donor_ax.title.set_text('Donor')
gut_donor_ax.set_ylabel("Secondary metabolite")

# gut fmt
att3 = samples_smMetadata.query('Type == "FMT"').droplevel('Type').T
sns.heatmap(att3,cmap='binary',ax=gut_patient_fmt_ax,cbar=False)
gut_patient_fmt_ax.title.set_text('FMT - Stool')

# gut placebo
att3 = samples_smMetadata.query('Type == "PLACEBO"').droplevel('Type').T
sns.heatmap(att3,cmap='binary',ax=gut_patient_placebo_ax, cbar_kws={'label': 'Mean relative abundance'})
gut_patient_placebo_ax.title.set_text('Placebo - Stool')

plt.tight_layout()
plt.savefig('results/antismash.pdf')
#plt.show()
