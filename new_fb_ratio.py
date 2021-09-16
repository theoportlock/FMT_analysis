#!/usr/bin/env python import numpy as np
import seaborn as sns
import numpy as np
import pandas as pd
import itertools
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

samples_metadata = pd.read_csv("newmergedmetadata.csv").set_index('ID_x')
gut_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
gut_msp_data = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
oral_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)
oral_msp_data = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)

taxa_type='phylum'
dtaxa_not_sum = gut_msp_data.join(gut_taxonomy[taxa_type], how='inner').set_index(taxa_type).T
otaxa_not_sum = oral_msp_data.join(oral_taxonomy[taxa_type], how='inner').set_index(taxa_type).T
otaxa = otaxa_not_sum.groupby(otaxa_not_sum.columns, axis=1).sum()
dtaxa = dtaxa_not_sum.groupby(dtaxa_not_sum.columns, axis=1).sum()

otaxa['Site'], dtaxa['Site'] = 'Oral', 'Gut'
joined = pd.concat([otaxa, dtaxa])
meta_joined = joined.join(samples_metadata[['Days after treatment', 'Type', 'id']], how='inner')
#meta_joined = meta_joined.groupby(meta_joined.index).first()
meta_joined['F/B'] = np.log2(meta_joined.Firmicutes / meta_joined.Bacteroidetes)
meta_joined = meta_joined.replace([np.inf, -np.inf], np.nan).dropna(subset=["F/B"], how="all")

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_palette(sns.color_palette("tab20").reverse())
fig, (gut_donor_ax, gut_patient_ax, oral_patient_ax) = plt.subplots(1, 3, sharey=True, gridspec_kw={'width_ratios': [1, 8, 8]})
x='Days after treatment'
y='F/B'
hue='Type'

donordf = meta_joined.query('Type == "DONOR"').drop('Days after treatment', axis=1).drop_duplicates()
donordf['F/B'] = np.log2(donordf.Firmicutes / donordf.Bacteroidetes)
sns.boxplot(data=donordf, y='F/B', palette=['C2'], linewidth=1, ax=gut_donor_ax)
sns.stripplot(data=donordf, y='F/B', palette=['C2'], jitter=False, linewidth=1, ax=gut_donor_ax)

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
plt.ylim(-0.5,12)
gut_donor_ax.title.set_text('Donor Stool')
gut_donor_ax.set_ylabel("$Log_2$(F/B ratio)")
gut_patient_ax.set_ylabel("")
gut_patient_ax.title.set_text('Stool')
oral_patient_ax.title.set_text('Saliva')
oral_patient_ax.set_ylabel("")
gut_patient_ax.legend([],[], frameon=False)
gut_donor_ax.legend([],[], frameon=False)

# Stats
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

#plt.show()
plt.tight_layout()
plt.savefig('results/FBratio.pdf')
