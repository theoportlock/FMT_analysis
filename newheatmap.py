#!/usr/bin/env python
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np

file_taxonomy="../downstream_data/taxo.csv"
file_meteor_data="../downstream_data/merged.final.mgs.med.vec.10M.csv"
file_names="../downstream_data/P15952_20_03_sample_info_stool.txt"
file_sample_type="../downstream_data/PROFITplaceboTab.csv"
file_day0="../downstream_data/PROFITbaseline.csv"
file_day7="../downstream_data/PROFITday7.csv"
file_day30="../downstream_data/PROFITday30.csv"
file_day90="../downstream_data/PROFITday90.csv"
file_siafile="../downstream_data/normalised.faecal.cytokines.markers.xlsx"

taxonomy = pd.read_csv(file_taxonomy, index_col=0)
data = pd.read_csv(file_meteor_data, index_col=0)
names = pd.read_table(file_names)
sample_type = pd.read_csv(file_sample_type, index_col=0)
day0=pd.read_csv(file_day0)
day7=pd.read_csv(file_day7)
day30=pd.read_csv(file_day30)
day90=pd.read_csv(file_day90)
siafile = pd.read_excel(file_siafile)

# add multi-index to data
sample_type = sample_type.set_index(sample_type.index.str.replace("P01","P").str.replace("\ .*", "", regex=True))
formatted_names = names['User ID'].str.split(expand=True)
formatted_names.loc[formatted_names[0]=="Donor", 2] = formatted_names.loc[formatted_names[0]=="Donor", 1]
formatted_names.loc[formatted_names[0]=="Donor", 0] = formatted_names.loc[formatted_names[0]=="Donor", 1]
formatted_names.iloc[:,2] = formatted_names.iloc[:,2].str.replace("Stool","0")
formatted_names['ID'] = names['NGI ID'].str.replace('P15952_','FMT')
formatted_names = pd.merge(formatted_names, sample_type, left_on=0, right_on="Patient", how='left')
formatted_names.fillna("DONOR", inplace=True)
data.columns = pd.MultiIndex.from_frame(formatted_names[[0,2,'Type','ID']], names=['Patient','Days after treatment', 'Type', 'ID'])
data.reset_index(inplace=True)
day0["Days after treatment"]=0
day7["Days after treatment"]=7
day30["Days after treatment"]=30
day90["Days after treatment"]=90
vals = pd.concat([day0,day7,day30,day90])
vals["Patient_Id"] = vals.Patient.str.replace("P01","P").str.replace("\ .*", "", regex=True)
siafile["id"] = siafile.id.str.replace("P01","P").str.replace("\ .*", "", regex=True)
siafile["time_point"] = siafile.time_point.str.replace("D","")
siafile['time_point'] = siafile.time_point.astype('int64')
merged_metadata = pd.merge(siafile, vals,  how='inner', left_on=['id','time_point'], right_on = ['Patient_Id','Days after treatment'])
#numeric_data = merged_metadata.select_dtypes(include=np.number)
fmerged_metadata = merged_metadata.drop(["Age at enrolment","MELD_y","Days after treatment"],axis=1)
ffmerged_metadata = fmerged_metadata.T[fmerged_metadata.isna().sum(axis=0) < 5].T
fgut_msp_data = gut_msp_data.T.reset_index()
fgut_msp_data["Days after treatment"] = fgut_msp_data["Days after treatment"].astype(int)
#normalized_df=(numeric_data-numeric_data.min())/(numeric_data.max()-numeric_data.min())
merged_MSP_metadata = pd.merge(fgut_msp_data, ffmerged_metadata, left_on=['Patient','Days after treatment'], right_on=['id','time_point'])
fmerged_MSP_metadata = merged_MSP_metadata.infer_objects()
numeric_data = fmerged_MSP_metadata.select_dtypes(include=np.number)
cor = numeric_data.corr(method='spearman')
fcor = cor.loc[cor.columns.str.contains('msp'),~cor.columns.str.contains('msp')]
fcor.dropna(axis=0,inplace = True)
ffcor = fcor
#ffcor = fcor[(fcor.min(axis=1) < -0.5) | (fcor.max(axis=1) > 0.5)]

#sns.heatmap(
#        ffcor,
#        cmap='coolwarm',
#        yticklabels=True,
#        xticklabels=True)

'''
get pearson pval
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

pval = numeric_data.corr(method=lambda x, y: pearsonr(x, y)[1]) - np.eye(len(cor.columns))
'''

sns.clustermap(
        ffcor,
                cmap='coolwarm',
                        yticklabels=False,
                                xticklabels=True)
