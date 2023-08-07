import pandas as pd
import numpy as np

edges = pd.read_csv('fsnewedges.csv')
nodes = pd.read_csv('nnnodes.csv')

nodes = nodes.loc[(nodes.VISIT == 'WEEK04') & (nodes.ARM == 'ACTIVE'),:].drop(['VISIT', 'ARM'], axis=1).T


edges['from'] = edges['from'].str.replace(r' \(.*', '', regex=True)
nodes.index = nodes.index.str.replace(r'  \(.*', '', regex=True)
elements = np.append(edges['from'].unique(),(edges['to'].unique()))
final = nodes.loc[nodes.index.isin(elements)]

final.to_csv('nnnnodes.csv')


edges.to_csv('fsnewedges.csv')
