## Import Statements
import numpy as np
import pandas as pd
import networkx as nx
import sys
import os
import OmicsIntegrator as oi
from math import log as math_log

## Anna Ritz 2020
## TODO add comments

def run(G,sources,targets,outprefix,dummy_edge_weight=5,edge_reliability=1,degree_penalty=3,prize=5,verbose=False):
	##TODO: this function writes a new interactome file and a prizes file, THEN runs PCSF. These should be moved outside this function.


	## OmicsIntegrator expects an undirected graph of 'ProteinA ProteinB Cost'
	interactome_file = prepare_interactome(G,outprefix)
	## defult params
	## defaults = {"w": 5, "b": 1, "g": 3, "edge_noise": 0.1, "dummy_mode": "terminals", "seed": 0, "skip_checks": False}
	## w: dummy edge weight
	## b: edge reliability
	## g: degree penalty
	graph = oi.Graph(interactome_file, {'w':dummy_edge_weight,'b':edge_reliability,'g':degree_penalty,'edge_noise':0.001, 'skip_checks':True})

	## OMicsIntegrator expets a prize file with 'Protein Prize' with a header row.
	prize_file = prepare_prizes(sources,targets,prize,outprefix)
	graph.prepare_prizes(prize_file)

	## run PCSF
	vertex_indices, edge_indices = graph.pcsf(verbosity_level=1)
	print('Forest has %d nodes and %d edges' %(len(vertex_indices),len(edge_indices)))

	## get forest
	forest, augmented_forest = graph.output_forest_as_networkx(vertex_indices, edge_indices)
	print('Objective value: %f' % (graph.pcsf_objective_value(forest)))

	return forest

def prepare_interactome(G,outprefix):
	print("preparing interactome....")
	outfile = '%s-interactome.txt' % (outprefix)
	## make undirected, taking smallest cost if edge is bidirected.
	seen = set()
	out = open(outfile,'w')
	out.write('protein1\tprotein2\tcost\n')
	for u,v in G.edges():
		if (u,v) in seen or (v,u) in seen:
			continue
		if G.has_edge(v,u):
			cost = min(G[u][v]['cost'],G[v][u]['cost'])
		else:
			cost = G[u][v]['cost']
		out.write('%s\t%s\t%f\n' % (u,v,cost))
		#out.write('%s\t%s\t%f\n' % (v,u,cost))
		seen.add((v,u))
		seen.add((u,v))
	out.close()
	print('wrote to %s' % (outfile))
	return outfile


def prepare_prizes(sources,targets,p,outprefix):
	print("preparing prizes...")
	outfile = '%s-prizes.txt' % (outprefix)
	out = open(outfile,'w')
	out.write('protein\tprize\n')
	for n in sources.union(targets):
		out.write('%s\t%f\n' % (n,p))
	out.close()
	print('wrote to %s' % (outfile))
	return outfile

## copied from PerfectLinker
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
            GR.add_edge(u,v,weight=float(w),cost=-math_log(max([0.000000001, w]))/math_log(10))
        else:
            u,v = e
            GR.add_edge(u,v)
    if verbose:
        print('length of edge set: {}'.format(len(set(edges))))
        print('number of edges in GR: {}'.format(len(list(GR.edges))))
    return GR

def get_labeled_nodes(fname: str):
    labels = pd.read_csv(fname,sep='\t')
    try:
        sources = list(labels[labels['Node type'] == 'source']['#Node'])
    except:
        sources = list(labels[labels['node_symbol'] == 'receptor']['#node'])
    try:
        sinks = list(labels[labels['Node type'] == 'target']['#Node'])
    except:
        sinks = list(labels[labels['node_symbol'] == 'tf']['#node'])
    return sources,sinks

def write_output(G,outfile,verbose=False):
	out = open(outfile,'w')
	out.write('#tail\thead\n')
	seen = set()
	for u,v in G.edges():
		if (v,u) in seen:
			continue
		seen.add((u,v))
		out.write('%s\t%s\n' % (u,v))
	out.close()
	if verbose:
		print('wrote to %s' % (outfile))
	return

def main(argv):
    interactome = argv[1]
    labeled_nodes = argv[2]
    dummy_edge_weight=float(argv[3])
    edge_reliability=float(argv[4])
    degree_penalty=float(argv[5])
    prize=float(argv[6])
    verbose = bool(argv[7])

    print('preprocessing inputs...')
    interactome = df_to_graph(interactome)
    sources,sinks = get_labeled_nodes(labeled_nodes)

    ## sources/sinks need to be sets here.
    sources = set(sources)
    sinks = set(sinks)

    outprefix = 'pcsf-w%d-b%d-g%d-p%d' % (int(dummy_edge_weight),int(edge_reliability),int(degree_penalty),int(prize))
    print('making prediction...')
    forest = run(interactome,sources,sinks,outprefix, \
    	dummy_edge_weight=dummy_edge_weight, \
    	edge_reliability=edge_reliability, \
    	degree_penalty=degree_penalty, \
    	prize=prize, verbose=verbose)

    print('saving prediction...')
    write_output(forest, '%s.csv' % (outprefix),verbose=verbose)

    return 

if __name__ == '__main__':
    main(sys.argv)



