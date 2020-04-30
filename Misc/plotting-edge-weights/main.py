#
# Tobias Rubel | rubelato@reed.edu
# Reed CompBio 
# 
# This program generates various plots of the edge weights 
# of interactions inside and outside of pathways of interest.
# ultimately, my goal is to create a node weighting scheme
# based on edge weights which makes HybridLinker better.
#
import os
import numpy as np
import pandas as pd


#first lets load the interactome
print('loading PPI...')
PPI = pd.read_csv('../../Interactomes/PathLinker_2018_human-ppi-weighted-cap0_75.txt',sep='\t')
#next lets get ahold of the edges in our pathways

PPATH = '../../Pathways'

pedges = []
for e in [x for x in os.listdir(PPATH) if 'edges' in x]:
    ne = pd.read_csv(os.path.join(PPATH,e),sep='\t').take([0,1],axis=1)
    pedges.append(ne)

pedges = pd.concat(pedges)
#filter out edges not in interactome and assign edge weights
Pathway_edges = pd.merge(PPI,pedges,on=['#tail','head'])

#get the complement too
keys = list(Pathway_edges.columns.values[:2])
i1 = PPI.set_index(keys).index
i2 = Pathway_edges.set_index(keys).index
NPathway = PPI[~i1.isin(i2)]

#compute statistics
pedge_weights = Pathway_edges['edge_weight'].mean()
print('mean weight of in pathway edges: {}'.format(pedge_weights))
PPI_edge_weights = PPI['edge_weight'].mean()
print('mean weight of arbitrary edges: {}'.format(PPI_edge_weights))
NPathway_edge_weights = NPathway['edge_weight'].mean()
print('mean weight of edges not in pathway: {}'.format(NPathway_edge_weights))

#now lets check score nodes by the edges they are part of.
#let the node weight = mean weight of edges containing node

access_node = lambda x: PPI[(PPI['#tail']==x) | (PPI['head']==x)]

weight_node = lambda x: access_node(x).mean()

print('building node dataframe...')
nlist = [x for y in [tuple(z) for z in PPI.take([0,1],axis=1).values] for x in y]
weight = []
i = 0
for x in nlist:
    i+=1
    weight.append(weight_node(x))
    print(i)

NDF = pd.DataFrame.from_dict({'node':nlist,'weight':weight})

NDF.to_csv('NDF.csv')
