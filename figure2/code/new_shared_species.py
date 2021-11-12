#!/usr/bin/env python
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

gmsp_samples = pd.read_csv("../../data/gutmsp.csv", index_col=0)
omsp_samples = pd.read_csv("../../data/oralmsp.csv", index_col=0).T

gmsp_gtaxonomy = pd.read_csv("../../../downstream_data/taxo.csv", index_col=-1)
omsp_otaxonomy = pd.read_csv("../../../oral_downstream_data/oraltaxo.csv", index_col=0)

taxaType='genus'
gtaxonomy_samples = gmsp_samples.join(gmsp_gtaxonomy[taxaType], how='inner').groupby(taxaType).sum() > 0
samples_gtaxonomy = gtaxonomy_samples.T

otaxonomy_samples = omsp_samples.join(omsp_otaxonomy[taxaType], how='inner').groupby(taxaType).sum() > 0
samples_otaxonomy = otaxonomy_samples.T


