#
# Tobias Rubel | rubelato@reed.edu
# CompBio
#
# This program is a control for neighbor based methods of pathway
# reconstruction. It is just DFS which leverages pathway information.
#
# It is a choice-point as to whether we demand that edges are part
# of paths which connect sources to sinks. I think there will have to 
# be two methods which do it.
#
#
import os
import sys

import pandas as pd
import numpy as np
import networkx as nx
from queue import LifoQueue

def df_to_graph(fname: str,verbose=True,weighted=True):
    """
    :fname   path/to/dataframe
    :returns nx graph
    """
    if weighted:
        df = pd.read_csv(fname,sep='\t').take([0,1,2],axis=1)
    else:
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
        if weighted:
            u,v,w = e
            GR.add_edge(u,v,weight=float(w))
        else:
            u,v = e
            GR.add_edge(u,v)
    if verbose:
        print('length of edge set: {}'.format(len(set(edges))))
        print('number of edges in GR: {}'.format(len(list(GR.edges))))
    return GR


def perfect_linker_bfs(v,G: nx.DiGraph,ground: nx.DiGraph,variety='nodes',verbose=False) -> nx.DiGraph:
    """
    :v          vertex
    :G          interactome
    :ground     ground truth network
    :returns    the largest subgraph of ground achievable through neighbor based methods

    :note       this method is guaranteed to have perfect precision if edges is chosen as variety.
    :fix?       it isn't totally clear that this is depth first
    """
    #initialize stack
    S = LifoQueue()
    S.put(v)
    #initialize discovered
    discovered = set()
    discovered.add(v)
    #initialize return graph
    ret = nx.DiGraph()
    ret.add_node(v)
    while not S.empty():
        if verbose:
            print('got {} out of {}'.format(len(discovered),len(ground.edges)))
        x = S.get()
        if verbose:
            print('processing node: {}'.format(x))
        f = lambda y: G[x][y]['weight']
        for n in sorted(list(G.neighbors(x)),key=f,reverse=True):
            if variety == 'nodes':
                if (not ((x,n) in discovered)) and (n in ground.nodes):
                    #add edge,vertex to ret
                    ret.add_node(n)
                    ret.add_edge(x,n)
                    discovered.add((x,n))
                    S.put(n)
            elif variety == 'edges':
                if (not (x,n) in discovered) and ((x,n) in ground.edges):
                    #add edge,vertex to ret
                    ret.add_node(n)
                    ret.add_edge(x,n)
                    discovered.add((x,n))
                    S.put(n)
    

    print('got {} out of {}'.format(len(discovered),len(ground.edges)))
    return ret

def src_snk(pathway: nx.DiGraph,labels: pd.DataFrame,verbose=True) -> nx.DiGraph:
    """
    :pathway graph
    :labels  dataframe of labels 
    :returns graph with source and sink added 
    """
    pathway.add_node('SRC')
    pathway.add_node('SNK')
    #get source nodes
    try:
        sources = list(labels[labels['Node type'] == 'source']['#Node'])
    except:
        sources = list(labels[labels['node_symbol'] == 'receptor']['#node'])
    if verbose:
        print('sources: {}'.format(sources))
    for s in sources:
        pathway.add_edge('SRC',s,weight=0)
    #get sink nodes
    try:
        sinks = list(labels[labels['Node type'] == 'target']['#Node'])
    except:
        sinks = list(labels[labels['node_symbol'] == 'tf']['#node'])
    if verbose:
        print('sinks: {}'.format(sinks))
    for s in sinks:
        pathway.add_edge(s,'SNK',weight=0)
    return pathway


def graph_to_file(G: nx.Graph,variety: str) -> None:
    df = pd.DataFrame({'#tail':[x[0] for x in G.edges],'head':[x[1] for x in G.edges],'KSP index':[x for x in range(1,len(G.edges)+1)]})
    df.to_csv('{}-PerfectLinker-DFS-Weighted.csv'.format(variety),index=False,sep='\t')


def main(argv):
    kind,interactome,pathway,labeled_nodes = argv[1:5]
    print('preprocessing inputs...')
    interactome = df_to_graph(interactome)
    pathway = df_to_graph(pathway,weighted=True)
    labeled_nodes = pd.read_csv(labeled_nodes,sep='\t')
    #add source and sink nodes to graph
    interactome = src_snk(interactome,labeled_nodes)
    pathway = src_snk(pathway,labeled_nodes)
    #print(pathway.edges)
    print('making prediction...')
    prediction = perfect_linker_bfs('SRC',interactome,pathway,kind)
    print('saving prediction...')
    graph_to_file(prediction,kind)

if __name__=="__main__":
    main(sys.argv)

    



                

             
 
