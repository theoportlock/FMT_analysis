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
numeric_data = merged_metadata.select_dtypes(include=np.number)
#normalized_df=(numeric_data-numeric_data.min())/(numeric_data.max()-numeric_data.min())

'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

cor_matrix = numeric_data.corr(method='spearman')
stocks = cor_matrix.index.values
cor_matrix = np.asmatrix(cor_matrix)

G = nx.from_numpy_matrix(cor_matrix)
G = nx.relabel_nodes(G,lambda x: stocks[x])
G.edges(data=True)

#function to create and display networks from the correlatin matrix.

def create_corr_network_1(G):
    #crates a list for edges and for the weights
    edges,weights = zip(*nx.get_edge_attributes(G,'weight').items())
    #positions
    positions=nx.circular_layout(G)
    #Figure size
    plt.figure(figsize=(15,15))
    #draws nodes
    nx.draw_networkx_nodes(G,positions,node_color='#DA70D6',
                           node_size=500,alpha=0.8)
    #Styling for labels
    nx.draw_networkx_labels(G, positions, font_size=8,
                            font_family='sans-serif')
    #draws the edges
    nx.draw_networkx_edges(G, positions, edgelist=edges,style='solid')
    # displays the graph without axis
    plt.axis('off')
    #saves image
    #plt.savefig("part1.png", format="PNG")
    plt.show()

create_corr_network_1(G)
'''
