#!/usr/bin/env python import numpy as np
import seaborn as sns
import numpy as np
import pandas as pd
import itertools
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

samples_metadata = pd.read_csv("newnewmetadata.csv").set_index('Sample ID')
samples_gutmsp = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0).T
samples_oralmsp = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0).T
samples_gutsum = pd.DataFrame(samples_gutmsp.agg(np.count_nonzero, axis=1), columns=['Species richness']).join(samples_metadata)

#samples_gutsum = pd.DataFrame(samples_gutmsp.agg(shannon, axis=1), columns=['Shannon Index']).join(samples_metadata)
samples_gutsum['Site'] = 'Gut'
samples_oralsum = pd.DataFrame(samples_oralmsp.agg(np.count_nonzero, axis=1), columns=['Species richness']).join(samples_metadata)
#samples_oralsum = pd.DataFrame(samples_oralmsp.agg(shannon, axis=1), columns=['Shannon Index']).join(samples_metadata)
samples_oralsum['Site'] = 'Oral'
meta_joined = pd.concat([samples_gutsum, samples_oralsum]).reset_index()
meta_joined['Type'] = meta_joined.Type.fillna('DONOR')
meta_joined.loc[meta_joined.Type == 'DONOR', 'Site'] = 'DONOR'

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_palette(sns.color_palette("tab20").reverse())
fig, (gut_donor_ax, gut_patient_ax, oral_patient_ax) = plt.subplots(1, 3, sharey=True, gridspec_kw={'width_ratios': [1, 8, 8]})
x='Days after treatment'
y='Species richness'
hue='Type'

donordf = meta_joined.query('Type == "DONOR"').drop('Days after treatment', axis=1).drop_duplicates()
sns.boxplot(data=donordf, y=y, palette=['C2'], linewidth=1, ax=gut_donor_ax)
sns.stripplot(data=donordf, y=y, palette=['C2'], jitter=False, linewidth=1, ax=gut_donor_ax)
#mannwhitneyu(donordf['Species richness'], meta_joined.loc[meta_joined['Days after treatment'] == 0]['Species richness'])
#mannwhitneyu(donordf['Species richness'], meta_joined.loc[(meta_joined['Days after treatment'] == 0) & (meta_joined['Type'] == 'FMT')]['Species richness'])

sns.boxplot(data=meta_joined.loc[meta_joined.Site == 'Gut'], x=x, y=y, hue=hue, linewidth=1, hue_order=['PLACEBO','FMT'], ax=gut_patient_ax)
sns.stripplot(data=meta_joined.loc[meta_joined.Site == 'Gut'], x=x, y=y, hue=hue, linewidth=1, jitter=False,hue_order=['PLACEBO','FMT'], ax=gut_patient_ax)
sns.pointplot(data=meta_joined.loc[(meta_joined.Site == 'Gut') & (meta_joined.Type == 'FMT')], x=x, y=y, hue='id', legend=False, scale=0.4, palette=['C1'], joined=True, ax=gut_patient_ax)
sns.pointplot(data=meta_joined.loc[(meta_joined.Site == 'Gut') & (meta_joined.Type == 'PLACEBO')], x=x, y=y, hue='id', legend=False, scale=0.4, palette=['C0'], joined=True, ax=gut_patient_ax)

sns.boxplot(data=meta_joined.loc[meta_joined.Site == 'Oral'], x=x, y=y, hue=hue, linewidth=1, hue_order=['PLACEBO','FMT'], ax=oral_patient_ax)
sns.stripplot(data=meta_joined.loc[meta_joined.Site == 'Oral'], x=x, y=y, hue=hue, linewidth=1, jitter=False, hue_order=['PLACEBO','FMT'], ax=oral_patient_ax)
sns.pointplot(data=meta_joined.loc[(meta_joined.Site == 'Oral') & (meta_joined.Type == 'FMT')], x=x, y=y, hue='id', legend=False, scale=0.4, palette=['C1'], joined=True, ax=oral_patient_ax)
sns.pointplot(data=meta_joined.loc[(meta_joined.Site == 'Oral') & (meta_joined.Type == 'PLACEBO')], x=x, y=y, hue='id', legend=False, scale=0.4, palette=['C0'], joined=True, ax=oral_patient_ax)

handles, labels = oral_patient_ax.get_legend_handles_labels()
plt.legend(handles[0:2], labels[0:2], title='Group', bbox_to_anchor=(1.001, 1), loc='upper left', fontsize='small')
#plt.ylim(-0.5,12)
gut_donor_ax.title.set_text('Donor Stool')
gut_donor_ax.set_ylabel("Species richness")
gut_patient_ax.set_ylabel("")
gut_patient_ax.title.set_text('Stool')
oral_patient_ax.title.set_text('Saliva')
oral_patient_ax.set_ylabel("")
gut_patient_ax.legend([],[], frameon=False)
gut_donor_ax.legend([],[], frameon=False)

# Stats
df = meta_joined.loc[(meta_joined.Site == 'Oral') & (meta_joined.Type != 'DONOR')]
box_pairs = list(itertools.combinations([(a,b) for a in df['Days after treatment'].unique() for b in df["Type"].unique()],2))
test_short_name = 'M.W.W'
pvalues = {}
for pair in box_pairs:
    data1 = df.groupby([x,hue])[y].get_group(pair[0])
    data2 = df.groupby([x,hue])[y].get_group(pair[1])
    stat, p = mannwhitneyu(data1, data2)
    pvalues[pair] = p
pvals = pd.Series(pvalues)
sigpvals = pvals.loc[pvals < 0.05]

ax, test_results = add_stat_annotation(
    oral_patient_ax,
    box_pairs=sigpvals.index,
    data=df,
    x=x,
    y=y,
    hue=hue,
    perform_stat_test=False,
    pvalues=sigpvals,
    test_short_name=test_short_name,
    text_format='full',
    #loc='outside',
    verbose=2)

df = meta_joined.loc[(meta_joined.Site == 'Gut') & (meta_joined.Type != 'DONOR')]
box_pairs = list(itertools.combinations([(a,b) for a in df['Days after treatment'].unique() for b in df["Type"].unique()],2))
test_short_name = 'M.W.W'
pvalues = {}
for pair in box_pairs:
    data1 = df.groupby([x,hue])[y].get_group(pair[0])
    data2 = df.groupby([x,hue])[y].get_group(pair[1])
    stat, p = mannwhitneyu(data1, data2)
    pvalues[pair] = p
pvals = pd.Series(pvalues)
sigpvals = pvals.loc[pvals < 0.05]

ax, test_results = add_stat_annotation(
    gut_patient_ax,
    box_pairs=sigpvals.index,
    data=df,
    x=x,
    y=y,
    hue=hue,
    perform_stat_test=False,
    pvalues=sigpvals,
    test_short_name=test_short_name,
    text_format='full',
    #loc='outside',
    verbose=2)

plt.show()
#plt.tight_layout()
#plt.savefig('results/newsigspec.pdf')
