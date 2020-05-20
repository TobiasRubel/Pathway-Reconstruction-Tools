#
# Tobias Rubel | rubelato@reed.edu
# CompBio
#
# this program implements PathLinker Plus, a modified PathLinker algorithm which
# uses two passes in order to get better performance.
# 
import sys
import os
import subprocess
import re
import pandas as pd
import networkx as nx
import random


#handle making new subgraph

def connect(vertices, G,verbose=False):
    """
    :vertices a set of vertices in V to connect
    :G        a graph to draw edges from
    :returns  a subgraph of G 
    """
    Q = nx.DiGraph(G.subgraph(list(vertices)))
    if verbose:
        print('number of edges in graph: {}'.format(len(G.edges)))
        print('number of edges in subgraph: {}'.format(len(Q.edges)))
    return nx.DiGraph(G.subgraph(list(vertices)))

def grow_neighbors(vertices,G):
    """
    :vertices a set of vertices in V
    :G        a graph to get neighbors from 
    :returns  a set of vertices plus their neighbors
    """
    return set(vertices).union({x for v in vertices for x in nx.all_neighbors(G,v)})


def gen_graph(vertices,G,k=1):
    print('number of vertices in new graph: {}'.format(len(vertices)))
    if k == 0:
        return connect(vertices,G)
    else:
        return gen_graph(grow_neighbors(vertices,G),G,k-1)

def gen_vertices(fname):
    df = pd.read_csv(fname,sep='\t').take([0,1],axis=1)
    return list(df.stack())

def df_to_graph(fname,verbose=False):
    """
    :fname   path/to/dataframe
    :returns nx graph
    """
    df = pd.read_csv(fname,sep='\t').take([0,1,2],axis=1)
    #df = pd.read_csv(fname,sep='\t',skiprows=0,names=['head','tail','weight']).take([0,1,2],axis=1)
    #df = df.drop(0)
    ff = df
    edges = [tuple(x) for x in df.values]
    #print(edges)
    vertices = list(df.stack())
    GR = nx.DiGraph()
    for v in vertices:
        GR.add_node(v)
    for e in edges:
        u,v,w = e
        GR.add_edge(u,v,weight=float(w))
    if verbose:
        print('length of edge set: {}'.format(len(set(edges))))
        print('number of edges in GR: {}'.format(len(list(GR.edges))))
    return GR

def graph_to_file(G: nx.Graph) -> None:
    df = pd.DataFrame({'#Node1':[x[0] for x in G.edges],'Node2':[x[1] for x in G.edges],'weight':[G[x[0]][x[1]]['weight'] for x in G.edges]})
    df.to_csv('tmp-interactome.csv',index=False,sep='\t')


def main(argv):
    """
    this routine should just pass through the arguments to pathlinker,
    then load the resulting nodes file and create a new interactome,
    then re-run pathlinker
    """
    #do first run
    os.chdir('Methods/PathLinker')
    os.system('rm tmp*')
    k,network,sources,pathway = argv[1],argv[2],argv[3],argv[4]
    #we need to modify network,pathway paths
    network = os.path.join('../../',network)
    pathway = os.path.join('../../',pathway)
    sources = os.path.join('../../',sources)
    argv = ['python2','run.py']+['-k',k]+['-o','tmp']+[network,sources,]
    print(argv)
    subprocess.call(argv)
    #generate new interactome
    pedges = next(x for x in os.listdir('.') if 'tmp' in x)
    V = gen_vertices(pedges)
    G = df_to_graph(network) 
    H = gen_graph(V,G,k=0)
    #print(H.edges)
    graph_to_file(H)
    #do second run
    print('Doing second run')
    os.chdir('../PerfectLinker-DFS-Weighted')
    argv_2 = ['python3','PL.py','nodes',network,'../PathLinker/tmp-interactome.csv',sources]
    subprocess.call(argv_2)
    os.replace('nodes-PerfectLinker-DFS-Weighted.csv','../../HybridLinker-DFS-Weighted.csv')

if __name__ == "__main__":
    main(sys.argv)


#testing stuff
"""
G = nx.DiGraph()
nodes=[1,2,3,4,5]
edges=[(1,2),(2,3),(3,4),(4,1),(4,5)]
for n in nodes:
    G.add_node(n)
for e in edges:
    G.add_edge(*e)

H = gen_graph({1,2,3},G,k=2)
"""

