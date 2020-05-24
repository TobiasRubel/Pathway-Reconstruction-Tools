import networkx as nx
import os
import pandas as pd
import sys
from math import log as math_log

## Computes two random walks with restarts: one from sources, and one from
## targets on reversed network. Computes the edge flux of each node and takes
## the product of the edge fluxes for each edge.
##TODO: add comments

## uses NetworkX's pagerank:
## https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.link_analysis.pagerank_alg.pagerank.html#networkx.algorithms.link_analysis.pagerank_alg.pagerank
def run(G, sources, targets, alpha=0.85, thres=0.80, verbose=False):

	## run RWR from sources
	print('forward flux calculation...')
	personalization = {s:1/len(sources) for s in sources}
	forward_pagerank = nx.pagerank(G,alpha=alpha,personalization=personalization,weight='weight')
	## We compute flux score for edge (u, v) by multiplying the visitation probability
	## of u by the edge weight and normalizing by the weighted out degree of u.
	for u in forward_pagerank:
		denom = len(list(G.successors(u)))
		for v in G.successors(u): # outgoing neighbors
			G[u][v]['forward_flux'] = forward_pagerank[u]*G[u][v]['weight']/denom

	print('backward flux calculation...')
	personalization = {t:1/len(targets) for t in targets}
	G_rev = G.reverse(copy=True)
	backward_pagerank = nx.pagerank(G_rev,alpha=alpha,personalization=personalization,weight='weight')
	## BACKWARD: We compute BACKWARDS flux score for edge (u, v) by multiplying the visitation probability
	## of v by the edge weight and normalizing by the weighted IN degree of v.
	for v in backward_pagerank:
		denom = len(list(G.predecessors(v)))
		for u in G.predecessors(v):
			G[u][v]['backward_flux'] = backward_pagerank[v]*G[u][v]['weight']/denom

	## combine flux
	print('combining fluxes...')
	tot_flux = 0
	for u,v in G.edges():
		#print(G[u][v]['forward_flux'],G[u][v]['backward_flux'])
		if G[u][v]['forward_flux'] == 0 or G[u][v]['backward_flux'] == 0:
			G[u][v]['neglog_combined_flux'] = sys.float_info.max ## maximum float value
			G[u][v]['combined_flux'] = 0
		else:
			G[u][v]['neglog_combined_flux'] = - (math_log(G[u][v]['forward_flux']) + math_log(G[u][v]['backward_flux']))
			G[u][v]['combined_flux'] = G[u][v]['forward_flux']*G[u][v]['backward_flux']
			tot_flux+=G[u][v]['combined_flux']

	print('sorting and taking predictions that comprise %f of total flux' % (thres))
	edges = {} #dictionary of (u,v): -log(flux*flux)
	running_sum = 0
	counter = 0
	## iterate through edges sorted by -log(flux*flux)
	for u,v,d in sorted(G.edges(data=True), key=lambda t: t[2]['neglog_combined_flux']):
		edges[(u,v)] = G[u][v]['neglog_combined_flux']
		running_sum+=G[u][v]['combined_flux']
		counter+=1
		if counter % 1000 == 0:
			print('%d of %d: %d of %d (%f)' % (counter,len(G.edges()),running_sum,tot_flux,running_sum/tot_flux))
		if running_sum/tot_flux > thres:
			print('theshold of %f limits predictions to %d edges'  %(thres,len(edges)))
			break

	return edges

## write the network
def write_output(edges,outfile,verbose=False):
	out = open(outfile,'w')
	out.write('#tail\thead\trank\n')
	for u,v in edges.keys():
		out.write('%s\t%s\t%f\n' % (u,v,edges[(u,v)]))
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

def main(argv):
	interactome = argv[1]
	labeled_nodes = argv[2]
	alpha = float(argv[3]) ## teleportation prob.
	thres = float(argv[4]) ## percentile threshold
	verbose = bool(argv[5])

	print('preprocessing inputs...')
	interactome = df_to_graph(interactome)
	sources,sinks = get_labeled_nodes(labeled_nodes)

	print('making prediction...')
	prediction = run(interactome, sources, sinks, alpha=alpha, thres=thres, verbose=verbose)

	print('saving prediction...')
	write_output(prediction,'RWR.csv',verbose=verbose)

	return

if __name__ == '__main__':
	main(sys.argv)
