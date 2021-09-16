#!/usr/bin/env python
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

samples_metadata = pd.read_csv("newmergedmetadata.csv").set_index('ID_x')
gut_msp_data = pd.read_csv("../downstream_data/merged.final.mgs.med.vec.10M.csv", index_col=0).T
oral_msp_data = pd.read_csv("../oral_merged_downstream_data/oralmsps.csv", index_col=0).T

variables=['Days after treatment', 'Type']
gut_meta_joined = gut_msp_data.join(samples_metadata[variables], how='inner')
oral_meta_joined = oral_msp_data.join(samples_metadata[variables], how='inner')

gut_meta_joined["Days after treatment"] = gut_meta_joined["Days after treatment"].fillna(0)
oral_meta_joined["Days after treatment"] = oral_meta_joined["Days after treatment"].fillna(0)

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set(font_scale=0.4)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
plt.rcParams["axes.grid"] = False

cmap='coolwarm'
var = ['Type', 'Days after treatment']

def plotter(df, v1, v2, axis):
    ndf = df.loc[(df["Days after treatment"] == v1) & (df['Type'] == v2)].drop(var, axis=1)
    sns.heatmap(ndf, ax=axis, cmap=cmap, vmax=0.00001, yticklabels=False, xticklabels=False, cbar=False, rasterized=True)
    axis.set_ylabel(f"{v2} D{v1}")

fig, (d, f0, f7, f30, f90, p0, p7, p30, p90) = plt.subplots(9, 2, constrained_layout=True)
d[1].axis('off')
d[1].set_title('Saliva')
d[0].set_title('Stool')
plotter(gut_meta_joined, 0, 'DONOR', d[0])
plotter(gut_meta_joined, 0, 'FMT', f0[0])
plotter(gut_meta_joined, 7, 'FMT', f7[0])
plotter(gut_meta_joined, 30, 'FMT', f30[0])
plotter(gut_meta_joined, 90, 'FMT', f90[0])
plotter(gut_meta_joined, 0, 'PLACEBO', p0[0])
plotter(gut_meta_joined, 7, 'PLACEBO', p7[0])
plotter(gut_meta_joined, 30, 'PLACEBO', p30[0])
plotter(gut_meta_joined, 90, 'PLACEBO', p90[0])
plotter(oral_meta_joined, 0, 'FMT', f0[1])
plotter(oral_meta_joined, 7, 'FMT', f7[1])
plotter(oral_meta_joined, 30, 'FMT', f30[1])
plotter(oral_meta_joined, 90, 'FMT', f90[1])
plotter(oral_meta_joined, 0, 'PLACEBO', p0[1])
plotter(oral_meta_joined, 7, 'PLACEBO', p7[1])
plotter(oral_meta_joined, 30, 'PLACEBO', p30[1])
plotter(oral_meta_joined, 90, 'PLACEBO', p90[1])

plt.tight_layout(pad=0.1)
#plt.show()
plt.savefig("results/barco.pdf", bbox_inches='tight', pad_inches=0.02)
