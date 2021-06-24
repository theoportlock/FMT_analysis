#!/usr/bin/env python
import numpy as np
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

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



baseline=pd.read_csv("../downstream_data/PROFITbaseline.csv",)
baseline["Days after treatment"]=0
day7=pd.read_csv("../downstream_data/PROFITday7.csv",)
day7["Days after treatment"]=7
day30=pd.read_csv("../downstream_data/PROFITday30.csv",)
day30["Days after treatment"]=30
day90=pd.read_csv("../downstream_data/PROFITday90.csv",)
day90["Days after treatment"]=90
sample_type = pd.read_csv("../downstream_data/PROFITplaceboTab.csv", index_col=0)
sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))

vals = pd.concat([baseline,day7,day30,day90])
vals["Patient_Id"] = vals.Patient.str.replace("P01","P").str.replace("\ .*", "", regex=True)
formatted_names = pd.merge(vals, sample_type, left_on="Patient_Id", right_on="Patient", how='left')
formatted_names = formatted_names[formatted_names.Type != np.NaN]
formatted_names = formatted_names.set_index(['Patient_Id','Days after treatment','Type']).reset_index()
formatted_names = formatted_names[formatted_names['Bilirubin'].notna()]
formatted_names['MELD_SCORE'] = formatted_names.apply(meld, axis=1)
formatted_names['CHILD'] = formatted_names.apply(child_pugh, axis=1)

fig, (meldax, childax) = plt.subplots(1, 2)
x='Days after treatment'
hue='Type'

sns.boxplot(data=formatted_names, x=x, y='MELD_SCORE', hue=hue, linewidth=1, ax=meldax)
sns.swarmplot(data=formatted_names, x=x, y='MELD_SCORE', hue=hue, linewidth=1, ax=meldax)
sns.boxplot(data=formatted_names, x=x, y='CHILD', hue=hue, linewidth=1, ax=childax)
sns.swarmplot(data=formatted_names, x=x, y='CHILD', hue=hue, linewidth=1, ax=childax)

#childax.legend([],[], frameon=False)
meldax.legend([],[], frameon=False)

sns.despine(trim=True, left=True)
box_pairs=[((0,'FMT'),(7,'FMT'))]

# Stats
ax, test_results = add_stat_annotation(
    meldax,
    data=formatted_names,
    x=x,
    y='MELD_SCORE',
    hue=hue,
    test='Mann-Whitney',
    text_format='full',
    #loc='outside',
    box_pairs=box_pairs,
    verbose=2)


plt.show()
plt.tight_layout()
#plt.savefig("results/meld_child.pdf")
