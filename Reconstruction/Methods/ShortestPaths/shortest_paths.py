import networkx as nx
import os
import pandas as pd
import sys

##TODO: add comments

## Computes the shortest paths between all receptor-TF pairs and returns
## the union of these paths as the reconstruction.
def run(G, sources, targets, verbose=False):
	edges = {} # dictionary of (u,v): number of paths
	for s in sources:
		## https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.shortest_paths.generic.shortest_path.html#networkx.algorithms.shortest_paths.generic.shortest_path
		shortest_path = nx.shortest_path(G,s,weight='cost')
		shortest_path_len = nx.shortest_path_length(G,s,weight='cost')
		for t in targets:
			path = shortest_path[t]
			for i in range(len(path)-1):
				u = path[i]
				v = path[i+1]
				if (u,v) not in edges:
					edges[(u,v)] = 1
				else:
					edges[(u,v)] += 1
		if verbose:
			print('source %s: stored %d edges' % (s,len(edges)))
	return edges

## write the network
def write_output(edges,outfile,verbose=False):
	out = open(outfile,'w')
	out.write('#tail\thead\tnum_paths\n')
	for u,v in edges.keys():
		out.write('%s\t%s\t%d\n' % (u,v,edges[(u,v)]))
	out.close()
	if verbose:
		print('wrote to %s' % (outfile))
	return

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
			GR.add_edge(u,v,weight=float(w))
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

def main(argv):
	interactome = argv[1]
	labeled_nodes = argv[2]
	verbose = bool(argv[3])

	print('preprocessing inputs...')
	interactome = df_to_graph(interactome)
	sources,sinks = get_labeled_nodes(labeled_nodes)

	print('making prediction...')
	prediction = run(interactome, sources, sinks, verbose=verbose)

	print('saving prediction...')
	write_output(prediction,'shortest_paths.csv',verbose=verbose)

	return

if __name__ == '__main__':
	main(sys.argv)
