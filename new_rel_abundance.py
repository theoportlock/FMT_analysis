#!/usr/bin/env python import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import random

samples_metadata = pd.read_csv("newmergedmetadata.csv").set_index('ID_x')
gut_taxonomy = pd.read_csv("../downstream_data/taxo.csv", index_col=0)
gut_msp_data = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0)
oral_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)
oral_msp_data = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)

# Format data
taxa_type='phylum'

# donor
dtaxa_not_sum = gut_msp_data.join(gut_taxonomy[taxa_type], how='inner').set_index(taxa_type).T
otaxa_not_sum = oral_msp_data.join(oral_taxonomy[taxa_type], how='inner').set_index(taxa_type).T
otaxa = otaxa_not_sum.groupby(otaxa_not_sum.columns, axis=1).sum()
dtaxa = dtaxa_not_sum.groupby(dtaxa_not_sum.columns, axis=1).sum()

otaxa['Site'], dtaxa['Site'] = 'Oral', 'Gut'
joined = pd.concat([otaxa, dtaxa])

meta_joined = joined.join(samples_metadata[['Days after treatment', 'Type']], how='inner')

sns.set(rc={'figure.figsize':(11.7,8.27)})
#sns.set_palette(sns.color_palette("Spectral", 18))
sns.set_palette(sns.color_palette("tab20").reverse())
fig, (gut_donor_ax, gut_patient_fmt_ax, gut_patient_placebo_ax, oral_patient_fmt_ax, oral_patient_placebo_ax) = plt.subplots(1, 5, sharey=True, gridspec_kw={'width_ratios': [1, 4, 4, 4, 4]})

donordf = meta_joined.query('Type == "DONOR"').drop('Days after treatment', axis=1).drop_duplicates()
donoraverage = donordf.set_index('Type').mean()
donornorm = pd.DataFrame(donoraverage/donoraverage.sum())
donornorm.T.plot.bar(stacked=True, ax=gut_donor_ax, legend=False)

gutfmtdf = meta_joined.query('Type == "FMT" & Site == "Gut"').drop('Type', axis=1)
gutfmtgrouped = gutfmtdf.groupby("Days after treatment").mean()
gutfmtnorm = gutfmtgrouped.apply(lambda x: x/x.sum(), axis=1)
gutfmtnorm.plot.bar(stacked=True, ax=gut_patient_fmt_ax, legend=False)

gutplacdf = meta_joined.query('Type == "PLACEBO" & Site == "Gut"').drop('Type', axis=1)
gutplacgrouped = gutplacdf.groupby("Days after treatment").mean()
gutplacnorm = gutplacgrouped.apply(lambda x: x/x.sum(), axis=1)
gutplacnorm.plot.bar(stacked=True, ax=gut_patient_placebo_ax, legend=False)

oralfmtdf = meta_joined.query('Type == "FMT" & Site == "Oral"').drop('Type', axis=1)
oralfmtgrouped = oralfmtdf.groupby("Days after treatment").mean()
oralfmtnorm = oralfmtgrouped.apply(lambda x: x/x.sum(), axis=1)
oralfmtnorm.plot.bar(stacked=True, ax=oral_patient_fmt_ax, legend=False)

oralplacdf = meta_joined.query('Type == "PLACEBO" & Site == "Oral"').drop('Type', axis=1)
oralplacgrouped = oralplacdf.groupby("Days after treatment").mean()
oralplacnorm = oralplacgrouped.apply(lambda x: x/x.sum(), axis=1)
oralplacnorm.plot.bar(stacked=True, ax=oral_patient_placebo_ax)

plt.legend(title=taxa_type, bbox_to_anchor=(1.001, 1), loc='upper left', fontsize='small')
plt.ylim(0,1)
gut_donor_ax.title.set_text('Donor')
gut_donor_ax.set_ylabel("Relative Abundance")
gut_patient_fmt_ax.title.set_text('FMT - Stool')
gut_patient_placebo_ax.title.set_text('Placebo - Stool')
oral_patient_fmt_ax.title.set_text('FMT - Saliva')
oral_patient_placebo_ax.title.set_text('Placebo - Saliva')

plt.show()
plt.tight_layout()
#plt.savefig('results/relative_abundance.pdf')
