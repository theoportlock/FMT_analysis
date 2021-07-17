#!/usr/bin/env python
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

%autoindent
file_taxonomy="../downstream_data/taxo.csv"
file_meteor_data="../downstream_data/merged.final.mgs.med.vec.10M.csv"
file_names="../downstream_data/P15952_20_03_sample_info_stool.txt"
file_sample_type="../downstream_data/PROFITplaceboTab.csv"
file_day0="../downstream_data/PROFITbaseline.csv"
file_day7="../downstream_data/PROFITday7.csv"
file_day30="../downstream_data/PROFITday30.csv"
file_day90="../downstream_data/PROFITday90.csv"
file_siafile="../downstream_data/normalised.faecal.cytokines.markers.xlsx"
file_gut_taxonomy="../downstream_data/taxo.csv"

taxonomy = pd.read_csv(file_taxonomy, index_col=0)
data = pd.read_csv(file_meteor_data, index_col=0)
names = pd.read_table(file_names)
sample_type = pd.read_csv(file_sample_type, index_col=0)
day0=pd.read_csv(file_day0)
day7=pd.read_csv(file_day7)
day30=pd.read_csv(file_day30)
day90=pd.read_csv(file_day90)
siafile = pd.read_excel(file_siafile)
gut_taxonomy = pd.read_csv(file_gut_taxonomy, index_col=0)

# add multi-index to data
sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))
formatted_names = names['User ID'].str.split(expand=True)
formatted_names.loc[formatted_names[0]=="Donor", 2] = formatted_names.loc[formatted_names[0]=="Donor", 1]
formatted_names.loc[formatted_names[0]=="Donor", 0] = formatted_names.loc[formatted_names[0]=="Donor", 1]
formatted_names.iloc[:,2] = formatted_names.iloc[:,2].str.replace("Stool","0")
formatted_names['ID'] = names['NGI ID'].str.replace('P15952_','FMT')
formatted_names = pd.merge(formatted_names, sample_type, left_on=0, right_on="Patient", how='left')
formatted_names.fillna("DONOR", inplace=True)

# Join taxo information
taxa_type='family'
taxa_not_sum = data.join(gut_taxonomy[taxa_type]).set_index(taxa_type).reset_index()
data = taxa_not_sum.groupby(taxa_type).mean().copy()
data.columns = pd.MultiIndex.from_frame(formatted_names[[0,2,'Type','ID']], names=['Patient','Days after treatment', 'Type', 'ID'])
data.reset_index(inplace=True)
data.drop('DONOR',level=2,axis='columns',inplace=True)

day0["Days after treatment"]=0
day7["Days after treatment"]=7
day30["Days after treatment"]=30
day90["Days after treatment"]=90
vals = pd.concat([day0,day7,day30,day90])
vals["Patient_Id"] = vals.Patient.str.replace("P01","P").str.replace("\ .*", "", regex=True)
siafile["id"] = siafile.id.str.replace("P01","P").str.replace("\ .*", "", regex=True)
siafile["time_point"] = siafile.time_point.str.replace("D","")
siafile['time_point'] = siafile.time_point.astype(int)
merged_metadata = pd.merge(siafile, vals,  how='inner', left_on=['id','time_point'], right_on = ['Patient_Id','Days after treatment'])
fmerged_metadata = merged_metadata.drop(["Age at enrolment","MELD_y","Days after treatment"],axis=1).infer_objects() 
ffmerged_metadata = fmerged_metadata.T[fmerged_metadata.isna().sum(axis=0) < 5].T.infer_objects()
fgut_msp_data = data.set_index(taxa_type).T.reset_index()
fgut_msp_data["Days after treatment"] = fgut_msp_data["Days after treatment"].astype(int)

merged_MSP_metadata = pd.merge(fgut_msp_data, ffmerged_metadata, left_on=['Patient','Days after treatment'], right_on=['id','time_point'])
fmerged_MSP_metadata = merged_MSP_metadata.infer_objects()
numeric_data = fmerged_MSP_metadata.dropna()._get_numeric_data()

def get_correlations(df):
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how="outer")
    correlations = dfcols.transpose().join(dfcols, how="outer")
    for ix, r in enumerate(df.columns):
        for jx, c in enumerate(df.columns):
            sp = spearmanr(df[r], df[c])
            correlations[c][r] = sp[0]
            pvalues[c][r] = sp[1] if ix > jx else np.nan  # Only store values below the diagonal
    return correlations.astype("float"), pvalues.astype("float")

correlations, uncorrected_p_values = get_correlations(numeric_data)
fcorrelations = correlations.loc[
    fgut_msp_data._get_numeric_data().columns[1:],
    ffmerged_metadata._get_numeric_data().columns]
funcorrected_p_values = uncorrected_p_values.loc[
    ffmerged_metadata.select_dtypes(include=np.number).columns,
    fgut_msp_data.select_dtypes(include=np.number).columns[1:]].T
ffuncorrected_p_values = funcorrected_p_values.loc[(fcorrelations.sum(axis=1) != 0), (fcorrelations.sum(axis=0) != 0)].copy()
ffcorrelations = fcorrelations.loc[(fcorrelations.sum(axis=1) != 0), (fcorrelations.sum(axis=0) != 0)].copy()
shape = ffuncorrected_p_values.values.shape
significant_matrix = multipletests(ffuncorrected_p_values.values.flatten())[0].reshape(shape)

g = sns.clustermap(
    ffcorrelations,
    cmap="vlag",
    vmin=-1,
    vmax=1, 
    yticklabels=True,
    xticklabels=True)

for tick in g.ax_heatmap.get_yticklabels():
    tick.set_rotation(0)

for i, ix in enumerate(g.dendrogram_row.reordered_ind):
    for j, jx in enumerate(g.dendrogram_col.reordered_ind):
        if i != j:
            text = g.ax_heatmap.text(
                j + 0.5,
                i + 0.5,
                "*" if significant_matrix[ix, jx] else "",
                ha="center",
                va="center",
                color="black",
            )
            text.set_fontsize(8)

plt.show()
plt.tight_layout()
#plt.savefig("clustermap.png", dpi=200)
