#!/usr/bin/env python
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation
from scipy.stats import mannwhitneyu
from skbio.diversity.alpha import shannon

sns.set(rc={'figure.figsize':(11.7,8.27)})

# Load data
samples_metadata = pd.read_csv('metadata.csv', index_col=0)
samples_gutmsp = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0).T
samples_oralmsp = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0).T

#samples_gutsum = pd.DataFrame(samples_gutmsp.agg(np.count_nonzero, axis=1), columns=['Species richness']).join(samples_metadata)
samples_gutsum = pd.DataFrame(samples_gutmsp.agg(shannon, axis=1), columns=['Shannon Index']).join(samples_metadata)
samples_gutsum['Site'] = 'Gut'
#samples_oralsum = pd.DataFrame(samples_oralmsp.agg(np.count_nonzero, axis=1), columns=['Species richness']).join(samples_metadata)
samples_oralsum = pd.DataFrame(samples_oralmsp.agg(shannon, axis=1), columns=['Shannon Index']).join(samples_metadata)
samples_oralsum['Site'] = 'Oral'
joined = pd.concat([samples_gutsum, samples_oralsum]).reset_index()
joined['Type'] = joined.Type.fillna('DONOR')
joined.loc[joined.Type == 'DONOR', 'Site'] = 'DONOR'

# Plot
x='Days after treatment'
y='Shannon Index'
hue='Type'
col='Site'

#g = sns.catplot(data=joined, x=x, y=y, col=col, hue=hue, linewidth=1, kind='box')
fig, (gut_donor_ax, gut_patient_ax, oral_patient_ax) = plt.subplots(1, 3, sharey=True, gridspec_kw={'width_ratios': [1, 7, 7]})
sns.set_theme(font_scale=0.8)
sns.boxplot(ax = gut_donor_ax, data=joined.loc[joined.Site == 'DONOR'], x="Type", y=y, linewidth=1).set_title("Donor stool")
sns.stripplot(ax = gut_donor_ax, data=joined.loc[joined.Site == 'DONOR'], x="Type", y=y, linewidth=1, jitter=False)
sns.boxplot(ax = gut_patient_ax, data=joined.loc[joined.Site == 'Gut'], x=x, y=y, hue=hue, linewidth=1).set_title("Patient stool")
sns.stripplot(ax = gut_patient_ax, data=joined.loc[joined.Site == 'Gut'], x=x, y=y, hue=hue, linewidth=1, jitter=False)
sns.pointplot(ax = gut_patient_ax, data=joined.loc[(joined.Site == 'Gut') & (joined.Type == 'FMT')], x=x, y=y, hue='Patient_ID', legend=False, scale=0.4, palette=['C1'], joined=True)
sns.pointplot(ax = gut_patient_ax, data=joined.loc[(joined.Site == 'Gut') & (joined.Type == 'PLACEBO')], x=x, y=y, hue='Patient_ID', legend=False, scale=0.4, palette=['C0'], joined=True)
sns.boxplot(ax = oral_patient_ax, data=joined.loc[joined.Site == 'Oral'], x=x, y=y, hue=hue, linewidth=1).set_title("Patient Saliva")
sns.stripplot(ax = oral_patient_ax, data=joined.loc[joined.Site == 'Oral'], x=x, y=y, hue=hue, linewidth=1, jitter=False)
sns.pointplot(ax = oral_patient_ax, data=joined.loc[(joined.Site == 'Oral') & (joined.Type == 'FMT')], x=x, y=y, hue='Patient_ID', legend=False, scale=0.4, palette=['C1'], joined=True)
sns.pointplot(ax = oral_patient_ax, data=joined.loc[(joined.Site == 'Oral') & (joined.Type == 'PLACEBO')], x=x, y=y, hue='Patient_ID', legend=False, scale=0.4, palette=['C0'], joined=True)

# Format Plot
gut_patient_ax.legend([],[], frameon=False)
gut_donor_ax.legend([],[], frameon=False)
gut_patient_ax.set_ylabel('')
oral_patient_ax.set_ylabel('')
gut_donor_ax.set_xlabel('')
sns.despine(trim=True, left=True)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# Stats
df = joined.loc[joined.Site == 'Gut']
box_pairs = itertools.combinations([(a,b) for a in df['Days after treatment'].unique() for b in df["Type"].unique()],2)
#box_pairs=[((0,'FMT'),(7,'FMT'))]
test_short_name = 'M.W.W'
pvalues = []
for pair in box_pairs:
    data1 = df.groupby([x,hue])[y].get_group(pair[0])
    data2 = df.groupby([x,hue])[y].get_group(pair[1])
    stat, p = mannwhitneyu(data1, data2)
    print(f"Performing {test_short_name} statistical test for equal variances on pair:",
          pair, "stat={:.2e} p-value={:.2e}".format(stat, p))
    pvalues.append(p)
print("pvalues:", pvalues)

ax, test_results = add_stat_annotation(
    gut_patient_ax,
    box_pairs=box_pairs,
    data=df,
    x=x,
    y=y,
    hue=hue,
    perform_stat_test=False,
    pvalues=pvalues,
    test_short_name=test_short_name,
    text_format='full',
    #loc='outside',
    verbose=2)

# Print
plt.tight_layout()
#plt.savefig("results/SI_species_richness.pdf", dpi=600)
plt.show()

'''
from scipy.stats import mannwhitneyu
df = joined.loc[(joined.Site=='Oral') & (joined.Type=='PLACEBO')]
#filt = pd.Series(gut.groupby('Patient_ID')['Days after treatment'].count() == 4)
filt = pd.Series(gut.groupby('Patient_ID')['Days after treatment'].count() == 3)
fdf = df[df.Patient_ID.isin(filt[filt].index)]
#oral = joined.loc[joined.Site=='Oral']
baseline = fdf.loc[fdf['Days after treatment'] == 0]
stats = fdf.groupby('Days after treatment').apply(lambda group: mannwhitneyu(baseline[y], group[y])[1])
statsdf = pd.DataFrame(stats.to_list(), columns=joined.columns, index=averageDays.index).drop(['Days after treatment'], axis=1)
'''
