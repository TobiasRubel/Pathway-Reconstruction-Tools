import sys
import time
import pandas as pd
import os

### GRAPHSPACE MODULES
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph

## networkk
import networkx as nx
from math import log as math_log

### UNIPROT MAPPING
import urllib.parse
import urllib.request

## https://stackoverflow.com/questions/35569042/ssl-certificate-verify-failed-with-python3
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

def connect_to_graphspace(username,password):
	graphspace = GraphSpace(username,password)
	return graphspace

def post(G,gs,gs_group):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
		if gs_group:
			gs.share_graph(graph=graph,group_name=gs_group)
	return graph

def post_to_graphspace(graphspace,G,edgefile,ground_truth_edges,sources,targets,title,gs_group,add_times,verbose,undirected=False,pos_only=False):
	edges = set()
	nodes = set()
	header = None

	ground_truth_nodes = {u for u,v in ground_truth_edges}.union({v for u,v in ground_truth_edges})

	edge_popup = {}
	with open(edgefile) as fin:
		for line in fin:
			if line[0] == '#':  # skip header
				header = line.strip().split()
				continue
			row = line.strip().split()
			u = row[0]
			v = row[1]
			if u == 'SRC' or v == 'SNK':
				continue
			if pos_only and not (u in ground_truth_nodes and v in ground_truth_nodes): # only take induced subgraph of postiive nodes
				continue
			edges.add((u,v))
			nodes.add(u)
			nodes.add(v)
			if u not in edge_popup:
				edge_popup[u] = {}
			if v not in edge_popup[u]:
				if undirected:
					edge_popup[u][v] = 'UNDIRECTED edge<br>'
				else:
					weight = G[u][v]['weight']
					cost = G[u][v]['cost']
					edge_popup[u][v] = '<b>Edge Weight</b>: %.4f<br><b>Edge Cost</b>: %.4f<br>' % (weight,cost)
				edge_popup[u][v] += '------<br>'
			for i in range(2,len(row)):
				edge_popup[u][v] += '<b>%s</b>: %s<br>' % (header[i],row[i])
			edge_popup[u][v] += '------<br>'
			if len(nodes) > 300 or len(edges) > 500:
				print('Truncating rankings!')
				title+= '(first %d nodes and %d edges)' % (len(nodes),len(edges))
				break

	uid2name = mapper(nodes)
	node_popup = {}
	for node in nodes:
		node_popup[node] = '<a href="https://www.uniprot.org/uniprot/%s" target="_blank">UniprotID: %s</a>' % (node,node)
		if node in sources:
			node_popup[node] += '<br>Receptor'
		if node in targets:
			node_popup[node] += '<br>Transcriptional Regulator'

	for u in edge_popup:
		for v in edge_popup[u]:
			edge_popup[u][v] = '%s -> %s <br>' % (uid2name.get(u,u),uid2name.get(v,v)) + edge_popup[u][v]

	G = GSGraph()
	if add_times:
		G.set_name('%s (TS %f)' % (title,time.time()))
	else:
		G.set_name(title)

	for node in nodes:

		if node in sources:
			color = '#609BBD'
			shape = 'diamond'
			height=80
		elif node in targets:
			color = '#F0BA75'
			shape = 'rectangle'
			valign='center'
			height=45
		else:
			color = '#C9C9C9'
			shape = 'ellipse'
			height=45
		if node in ground_truth_nodes:
			#print(node,'IN GROUND TRUTH')
			border_width=4
			border_color='#000000'
		else:
			border_width=2
			border_color = '#ACACAC'
		G.add_node(node,label=uid2name.get(node,node),popup=node_popup[node])
		G.add_node_style(node,shape=shape,color=color,width=90,height=height,border_width=border_width,border_color=border_color)

	added_edges = set()
	for u,v in edges:
		if (u,v) in ground_truth_edges or (v,u) in ground_truth_edges:
			width=4
			color='#000000'
		else: # ignored or negatives
			width=2
			color='#ACACAC'
		if undirected:
			G.add_edge(u,v,popup=edge_popup[u][v])
			G.add_edge_style(u,v,width=width,color=color)
		else:
			if (v,u) in edges:
				if (v,u) in added_edges: # already added - skip.
					continue
				else: # add bidirected edge. Right now this looks like an undirected edge.
					## get other popup too
					bidirected_popup = 'Bidirected Edge<br>----------------<br><br>' + edge_popup[u][v] + '<br>' + edge_popup[v][u]
					G.add_edge(u,v,popup=bidirected_popup)
					G.add_edge_style(u,v,width=width,color=color)
					added_edges.add((u,v))
			else: ## uni-directed edge. Add.
				G.add_edge(u,v,directed=True,popup=edge_popup[u][v])
				G.add_edge_style(u,v,directed=True,width=width,color=color)
	graph = post(G,graphspace,gs_group)
	if verbose:
		print('posted graph "%s" and shared with group "%s"' % (title,gs_group))

	return

def mapper(nodes):
	url = 'https://www.uniprot.org/uploadlists/'

	params = {
	'from': 'ACC+ID',
	'to': 'GENENAME',
	'format': 'tab',
	'query': ' '.join(list(nodes))
	}

	data = urllib.parse.urlencode(params)
	data = data.encode('utf-8')
	req = urllib.request.Request(url, data)
	with urllib.request.urlopen(req) as f:
	   response = f.read()
	answer = response.decode('utf-8')

	uid2name = {}
	for line in answer.split('\n'):
		if 'From' in line: # skip header
			continue
		row = line.strip().split()
		if len(row)==2:
			uid2name[row[0]] = row[1]
	return uid2name

## copied from PerfectLinker
def df_to_graph(fname: str,verbose=True,weighted=True):
	"""
	:fname   path/to/dataframe
	:returns nx graph
	"""
	print('interactome:',fname)
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

def load_df_tab(name:str):
    return pd.read_csv(name,sep='\t',engine='python')

def main(args):
	print('ARGUMENTS:',args[1:])
	username,password,name,directory,draft = args[1],args[2],args[3],args[4],eval(args[5])
	print('Draft?',draft)
	if draft:
		gs_group = None
		add_times = True
		verbose = True
	else:
		gs_group = 'reconstruction-traversals'
		add_times = False
		verbose=False
	## connect to GS server
	graphspace = connect_to_graphspace(username,password)

	## get relevant inputs
	print('preprocessing interactome...')
	#interactome = load_df_tab(os.path.join(directory,'interactome.csv'))
	interactome = df_to_graph(os.path.join(directory,'interactome.csv'))

	print('get ground truth...')
	ground = load_df_tab(os.path.join(directory,'ground.csv'))
	ground_truth_edges = {(x[0],x[1]) for x in ground.values if x[0] != 'SRC' and x[1] != 'SNK'}
	sources,sinks = get_labeled_nodes(os.path.join(directory,'ground-nodes.csv'))

	print('set prediction file...')
	predfile = os.path.join(directory,'ranked-edges.csv')

	post_to_graphspace(graphspace,interactome,predfile,ground_truth_edges,sources,sinks,name,gs_group,add_times,verbose)
	post_to_graphspace(graphspace,interactome,predfile,ground_truth_edges,sources,sinks,name + ' (ground truth nodes)',gs_group,add_times,verbose,pos_only=True)

if __name__ == '__main__':
	main(sys.argv)
