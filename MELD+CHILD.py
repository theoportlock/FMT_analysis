#!/usr/bin/env python
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

#%autoindent
#3.78×ln[serum bilirubin (mg/dL)] + 11.2×ln[INR] + 9.57×ln[serum creatinine (mg/dL)] + 6.43k
def meld(df):
    from numpy import log as ln
    return 3.78 * ln(df['Bilirubin']*10**-6) + 11.2 * ln(df['INR']) + 9.57 * ln(df['Creatinine']) + 6.43

def child_pugh(df):
    result = 0
    if df["Bilirubin"] > 34:
        if df["Bilirubin"] > 50:
            result += 3
        else:
            result += 2
    else:
        result += 1
    if df["Albumin"] / 10 < 3.5:
        if df["Albumin"] / 10 < 2.8:
            result += 3
        else:
            result += 2
    else:
        result += 1
    if df['INR'] > 1.7:
        if df['INR'] > 2.3:
            result += 3
        else:
            result += 2
    else:
        result += 1
    if df['Ascites'] > 0:
        if df['Ascites'] > 1:
            result += 3
        else:
            result += 2
    else:
        result += 1
    if df['HE grade'] > 0:
        if df['HE grade'] > 1:
            result += 3
        else:
            result += 2
    else:
        result += 1
    return result


formatted_names = pd.read_csv("metadata.csv", index_col=0)
formatted_names['MELD_SCORE'] = formatted_names.apply(meld, axis=1)
formatted_names['CHILD'] = formatted_names.apply(child_pugh, axis=1)

sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_theme(font_scale=0.8)

fig, (meldax, childax) = plt.subplots(1, 2)
x='Days after treatment'
y1='MELD_SCORE'
y2='CHILD'
hue='Type'

sns.boxplot(data=formatted_names, x=x, y=y1, hue=hue, linewidth=1, ax=meldax).set_title("Patient stool")
sns.stripplot(data=formatted_names, x=x, y=y1, hue=hue, linewidth=1, jitter=False, ax=meldax)
sns.pointplot(data=formatted_names.loc[formatted_names.Type == 'FMT'], x=x, y=y1, hue='Patient_ID', legend=False, scale=0.4, palette=['C1'], joined=True, ax=meldax)
sns.pointplot(data=formatted_names.loc[formatted_names.Type == 'PLACEBO'], x=x, y=y1, hue='Patient_ID', legend=False, scale=0.4, palette=['C0'], joined=True, ax=meldax)

sns.boxplot(data=formatted_names, x=x, y=y2, hue=hue, linewidth=1, ax=childax).set_title("Patient stool")
sns.stripplot(data=formatted_names, x=x, y=y2, hue=hue, linewidth=1, jitter=False, ax=childax)
sns.pointplot(data=formatted_names.loc[formatted_names.Type == 'FMT'], x=x, y=y2, hue='Patient_ID', legend=False, scale=0.4, palette=['C1'], joined=True, ax=childax)
sns.pointplot(data=formatted_names.loc[formatted_names.Type == 'PLACEBO'], x=x, y=y2, hue='Patient_ID', legend=False, scale=0.4, palette=['C0'], joined=True, ax=childax)

#childax.legend([],[], frameon=False)
meldax.legend([],[], frameon=False)
plt.legend(bbox_to_anchor=(1.001, 1), loc='upper left', fontsize='small')

sns.despine(trim=True, left=True)

#plt.show()
plt.tight_layout()
plt.savefig("results/meld_child.pdf")
