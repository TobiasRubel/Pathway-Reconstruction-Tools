#
# Tobias Rubel | rubelato@reed.edu
# Reed CompBio
#
# This program correlates number of sources and sinks to the number of edges/nodes
#
import os
import pandas as pd

DATA_PATH = '/home/tobias/Documents/Work/CompBio/ritz-data/pathways'

def load_df(name):
    global DATA_PATH
    return pd.read_csv(os.path.join(DATA_PATH,name),sep='\t')


def src_nodes(node_df):
    sources = list(node_df[node_df['node_symbol'] == 'receptor']['#node'])
    return sources 

def snk_nodes(node_df):
    sinks = list(node_df[node_df['node_symbol'] == 'tf']['#node'])
    return sinks

def len_pathway(edge_df):
    return len(edge_df['#tail'])

def len_pathway_nodes(edge_df):
    nodes = set(edge_df['#tail']).union(set(edge_df['head']))
    return len(nodes)


def len_labeled_nodes(node_df):
    return len(src_nodes(node_df))+len(snk_nodes(node_df))

def fetch_arguments():
    global DATA_PATH
    pathways = [set((x,y)) for x in os.listdir(DATA_PATH) for y in os.listdir(DATA_PATH) if x.split('-')[:-1] == y.split('-')[:-1] and x != y]
    sorted_pathways = [sorted(tuple(x),key = lambda x: x.split('-')[-1]) for x in pathways]
    processed_pathways = [(os.path.join(DATA_PATH,x),os.path.join(DATA_PATH,y)) for (x,y) in sorted_pathways]
    return list(set(processed_pathways))


def main():
    ARGS = fetch_arguments()
    l, e = [],[]
    for arg in ARGS:
        edges,nodes = arg
        l.append(len_labeled_nodes(load_df(nodes)))
        e.append(len_pathway(load_df(edges)))
    df = pd.DataFrame({'labels':l,'edges':e})
    df.to_csv('edge_vs_src+snk.csv')


if __name__ == "__main__":
    main()
