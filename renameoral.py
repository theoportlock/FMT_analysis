file_oral_msp_data="../oral_merged_downstream_data/merged.final.mgs.med.vec.10M.csv"
sm = pd.read_csv('metadata.csv')
oral_formatted_names = oral_names['User ID'].str.split(expand=True)
oral_formatted_names.loc[oral_formatted_names[2]=="raw", 2] = 0
oral_formatted_names[2] = oral_formatted_names[2].astype(int)
oral_formatted_names.merge(sm, right_on=['Patient_ID', 'Days after treatment'], left_on=[0,2], how='left')
ofn = oral_formatted_names.merge(sm, right_on=['Patient_ID', 'Days after treatment'], left_on=[0,2], how='left')
oral_msp_data.columns = pd.MultiIndex.from_frame(
        ofn[[0,2,'Type','ID']],
                names=['Patient_ID','Days after treatment', 'Type', 'ID'])
oral_msp_data
oral_msp_data.T
nomd = oral_msp_data.T.reset_index()
fnomd = nomd.iloc[:,3:]
fnomd
ffnomd = fnomd.set_index('ID')
fffnomd = ffnomd.loc[ffnomd.index.dropna()].T
fffnomd.to_csv('oralmsps.csv')
