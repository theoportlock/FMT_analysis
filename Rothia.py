#!/usr/bin/env python
import pandas as pd
import seaborn as sns

samples_metadata = pd.read_csv('metadata.csv', index_col=0)
msp_samples = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0)
msp_taxonomy = pd.read_csv("../oral_downstream_data/oraltaxo.csv", index_col=0)

taxaType='species'
taxonomy_samples = msp_samples.join(msp_taxonomy[taxaType], how='inner').groupby(taxaType).sum()
samples_taxonomy = taxonomy_samples.T
samples_taxonomyMetadata = samples_taxonomy.join(samples_metadata, how='inner')

tax = 'Centipeda peridontii'
tax = 'Rothia dentocariosa'
tax = 'Streptococcus gordonii'
sns.lmplot(data = samples_taxonomyMetadata, x=tax, y='MELD', hue='Type')
plt.show()
