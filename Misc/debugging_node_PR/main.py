#
# Tobias Rubel | rubelato@reed.edu
# Reed CompBio
# 
# This program is trying to figure out why
# HybridLinker is having weird node PR results.
# presumably the top-k value of HybridLinker
# should always align perfectly with the top-k
# value of PathLinker...
#

import pandas as pd

#first set which pathway we are looking at
print('loading HL and PL...')
HL = pd.read_csv('../../Validation/PR/data/HybridLinker_PathLinker_2018_human-ppi-weighted-cap0_75_TGF_beta_Receptor_500/ranked-edges.csv',sep='\t')

PL = pd.read_csv('../../Validation/PR/data/PathLinker_PathLinker_2018_human-ppi-weighted-cap0_75_TGF_beta_Receptor_500/ranked-edges.csv',sep='\t')

#print('HL:\n{}'.format(HL))

#print('PL:\n{}'.format(PL))

#HybridLinker shouldn't have any nodes not in PathLinker:
print('fetching their respective node sets...')
nodes = lambda df: {x for x in df.take([0,1],axis=1).stack().values}

HL_Nodes = nodes(HL)

PL_Nodes = nodes(PL)

#print('HL Nodes:\n{}'.format(HL_Nodes))

#print('PL Nodes:\n{}'.format(PL_Nodes))

Wrong_Nodes = {x for x in HL_Nodes if not x in PL_Nodes}

print('There should be no nodes in HL that are not in PL.\nthe following nodes were found. {}'.format(Wrong_Nodes))
