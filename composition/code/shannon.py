#!/usr/bin/env python
import numpy as np
import seaborn as sns
import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statannot import add_stat_annotation
from skbio.diversity.alpha import shannon

samples_metadata = pd.read_csv("../../data/newmergedmetadata.csv").set_index('ID_x')
samples_gutmsp = pd.read_csv("../../data/gutmsp.csv", index_col=0).T
samples_oralmsp = pd.read_csv("../../data/oralmsp.csv", index_col=0).T

samples_gutsum = pd.DataFrame(samples_gutmsp.agg(shannon, axis=1), columns=['Shannon Diversity Index']).join(samples_metadata)
samples_gutsum['Site'] = 'Gut'
samples_oralsum = pd.DataFrame(samples_oralmsp.agg(shannon, axis=1), columns=['Shannon Diversity Index']).join(samples_metadata)
samples_oralsum['Site'] = 'Oral'

meta_joined = pd.concat([samples_gutsum, samples_oralsum]).reset_index()
meta_joined['Type'] = meta_joined.Type.fillna('DONOR')
meta_joined.loc[meta_joined.Type == 'DONOR', 'Site'] = 'DONOR'

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_palette(sns.color_palette("tab20").reverse())
sns.set_style("whitegrid")
fig, (gut_donor_ax, gut_patient_ax, oral_patient_ax) = plt.subplots(1, 3, sharey=True, gridspec_kw={'width_ratios': [1, 8, 8]})
x='Days after treatment'
y='Shannon Diversity Index'
hue='Type'

donordf = meta_joined.query('Type == "DONOR"').drop('Days after treatment', axis=1).drop_duplicates()

sns.boxplot(data=donordf, y=y, palette=['C2'], linewidth=1, ax=gut_donor_ax)
sns.boxplot(data=meta_joined.loc[meta_joined.Site == 'Gut'], x=x, y=y, hue=hue, linewidth=1, hue_order=['PLACEBO','FMT'], ax=gut_patient_ax)
sns.boxplot(data=meta_joined.loc[meta_joined.Site == 'Oral'], x=x, y=y, hue=hue, linewidth=1, hue_order=['PLACEBO','FMT'], ax=oral_patient_ax)

handles, labels = oral_patient_ax.get_legend_handles_labels()
plt.legend(handles[0:2], labels[0:2], title='Group', bbox_to_anchor=(1.001, 1), loc='upper left', fontsize='small')
gut_donor_ax.title.set_text('Donor Stool')
gut_donor_ax.set_ylabel(y)
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

try:
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
        text_format='star',
        verbose=2)
except:
    pass

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

try:
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
        text_format='star',
        verbose=2)
except:
    pass

plt.tight_layout()
plt.savefig('../results/shannon.pdf')
#plt.show()
