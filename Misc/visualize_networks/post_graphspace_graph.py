#
# Tobias Rubel | rubelato@reed.edu
# CompBio
#
# This program is designed to visualize networks and predictions of them using graphspace. 
# it can be run via python3 main.py [network-name]
#
import os
import sys
import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns
#import primefac
from itertools import chain, combinations
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph

#declare global variables
graphspace = GraphSpace('rubelato@reed.edu','password')

#where all the pathways reside 
PATHWAY_PATH = '/home/tobias/Documents/Work/CompBio/ritz-data/pathways'




def load_df_tab(name:str):
    return pd.read_csv(name,sep='\t')

def df_to_graph(fname,verbose=False):
    """
    :fname   path/to/dataframe
    :returns nx graph
    """
    df = pd.read_csv(fname,sep='\t').take([0,1],axis=1)
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
        u,v = e
        GR.add_edge(u,v)
    if verbose:
        print('length of edge set: {}'.format(len(set(edges))))
        print('number of edges in GR: {}'.format(len(list(GR.edges))))
    return GR

def post_network(G,colo,v_colors,e_colors,color_key):
    """
    :G Graph 1
    :H Graph 2
    """
    GS = GSGraph()
    GS.set_name('TWEAK-Ground-vs-PerfectLinker-edges')
    GS.set_tags(['CompBio'])
    for v in G.nodes:
        GS.add_node(v,label=v)
        GS.add_node_style(v,height=50,width=50,color=colo[v_colors[v]])
    for e in G.edges:
        GS.add_edge(*e)
        GS.add_edge_style(*e,directed=True,color=colo[e_colors[e]])
    GS.set_data(data={'description': '{}'.format(color_key)})
    post(GS,graphspace)

def gen_colors(G,H):
    pass


def is_prime(n):
    """"pre-condition: n is a nonnegative integer
    post-condition: return True if n is prime and False otherwise."""
    if n < 2: 
         return False;
    if n % 2 == 0:             
         return n == 2  # return False
    k = 3
    while k*k <= n:
         if n % k == 0:
             return False
         k += 2
    return True

def primes(n):
    lat = []
    i = 2
    while len(lat) < n:
        if is_prime(i):
            lat.append(i)
        i += 1
    return lat

def gunion(lat):
    G = nx.DiGraph()
    for H in lat:
        V = H.nodes
        E = H.edges
        for v in V:
            try:
                G.add_node(v)
            except:
                print('node {} already in G'.format(v))
        for e in E:
            try:
                G.add_edge(*e)
            except:
                print('edge {} already in G'.format(e))
    return G

def prep_network(lat: list,labels:list):
    """
    :lat list of nx graphs
    """
    #we use these ids in order to index what colors to map things to
    ids = dict(zip(lat,primes(len(lat))))
    vertices = dict()
    for G in lat:
        for n in G.nodes():
            try:
                vertices[n] = vertices[n]*ids[G]
            except:
                vertices[n] = ids[G]
    color_map = colorize(ids)
    edges = dict()
    for G in lat:
        for n in G.edges():
            try:
                edges[n] = edges[n]*ids[G]
            except:
                edges[n] = ids[G]
    G = gunion(lat)
    return G,color_map,vertices,edges,dict(zip(powerset(labels),color_map.values()))

def product(list):
    p = 1
    for i in list:
        p *= i
    return p

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def colorize(idict):
    """
    :idict dict of graphs -> primes which serves as a key 
    """
    ncolors = 2**len(idict)
    palette = sns.color_palette(n_colors=ncolors).as_hex()
    all_ids = [product(i) for i in powerset(list(idict.values()))]
    id_to_color = dict(zip(all_ids,palette))
    return id_to_color


def post(G, gs):
    try:
        graph = gs.update_graph(G)
    except:
        graph = gs.post_graph(G)
    return graph

def main(argv):
    #G = df_to_graph(argv[1])
    #H = df_to_graph(argv[2])
    graphs = [df_to_graph(G) for G in argv[1:]]
    M = prep_network(graphs,argv[1:])
    post_network(*M)


if __name__ == "__main__":
    main(sys.argv)


