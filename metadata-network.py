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
siafile = pd.read_excel(file_siafile, engine = 'openpyxl')
othermeta = pd.read_csv('../downstream_data/Donorandpatientbaselinecharacteristics.tsv', sep='\t')

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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

cor_matrix = numeric_data.corr(method='spearman')
edges = cor_matrix.stack().reset_index()
edges.columns = ['asset_1','asset_2','correlation']
edges = edges.loc[edges['asset_1'] != edges['asset_2']].copy()

#create undirected graph with weights corresponding to the correlation magnitude
G0 = nx.from_pandas_edgelist(edges, 'asset_1', 'asset_2', edge_attr=['correlation'])

#print out the graph info
#check number of nodes and degrees are as expected (all should have degree = 38, i.e. average degree = 38)
print(nx.info(G0))

# 'winner takes all' method - set minium correlation threshold to remove some edges from the diagram
threshold = 0.5

# create a new graph from edge list
Gx = nx.from_pandas_edgelist(edges, 'asset_1', 'asset_2', edge_attr=['correlation'])

# list to store edges to remove
remove = []
# loop through edges in Gx and find correlations which are below the threshold
for asset_1, asset_2 in Gx.edges():
    corr = Gx[asset_1][asset_2]['correlation']
    #add to remove node list if abs(corr) < threshold
    if abs(corr) < threshold:
        remove.append((asset_1, asset_2))

# remove edges contained in the remove list
Gx.remove_edges_from(remove)
print(str(len(remove)) + " edges removed")

def assign_colour(correlation):
    if correlation <= 0:
        return "#ffa09b"  # red
    else:
        return "#9eccb7"  # green

def assign_thickness(correlation, benchmark_thickness=2, scaling_factor=3):
    return benchmark_thickness * abs(correlation)**scaling_factor


def assign_node_size(degree, scaling_factor=50):
    return degree * scaling_factor


# assign colours to edges depending on positive or negative correlation
# assign edge thickness depending on magnitude of correlation
edge_colours = []
edge_width = []
for key, value in nx.get_edge_attributes(Gx, 'correlation').items():
    edge_colours.append(assign_colour(value))
    edge_width.append(assign_thickness(value))

# assign node size depending on number of connections (degree)
node_size = []
for key, value in dict(Gx.degree).items():
    node_size.append(assign_node_size(value))

#create minimum spanning tree layout from Gx (after small correlations have been removed)
mst = nx.minimum_spanning_tree(Gx)

edge_colours = []

#assign edge colours
for key, value in nx.get_edge_attributes(mst, 'correlation').items():
    edge_colours.append(assign_colour(value))


#draw minimum spanning tree. Set node size and width to constant
nx.draw(mst, with_labels=True, pos=nx.fruchterman_reingold_layout(mst),
        node_size=200, node_color="#e1575c", edge_color=edge_colours,
       width = 1.2)

#set title
plt.title("Asset price correlations - Minimum Spanning Tree",fontdict=font_dict)
plt.show()



'''
formatted_names[2] = formatted_names[2].astype(int)
newmetadata = pd.merge(merged_metadata, formatted_names, left_on=['id', 'time_point'], right_on=[0,2], how='outer')
othermeta = pd.read_csv('../downstream_data/Donorandpatientbaselinecharacteristics.tsv', sep='\t')
othermeta['ID'] = othermeta.ID.str.replace('P01','P')
othermeta['ID'] = othermeta.ID.str.replace(' .*','')
othermeta['ID'] = othermeta.ID.str.replace('B','').str.replace('A','').str.replace('C','')
othermeta['ID'] = othermeta.ID.str.replace('D','')
othermeta = othermeta.groupby('ID').first().reset_index()
mmdata = pd.merge(newmetadata, othermeta, left_on=['id'], right_on=['ID'], how='outer')
mmdata.to_csv('newmergedmetadata.csv')
'''
