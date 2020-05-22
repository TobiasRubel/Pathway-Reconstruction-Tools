## makes latex-formatted table of pathways.
import glob
import pandas as pd

PATHWAY_DIR = '../../Pathways/'

FULL_NAMES = {'TSLP':'Thymic stromal lymphopoietin',
'ID':'Inhibitor of differentiation',
'KitReceptor':'Kit Receptor',
'IL6':'Interleukin-6',
'AndrogenReceptor':'Androgen Receptor',
'FSH':'Follicle-stimulating hormone',
'Prolactin':'Prolactin',
'Hedgehog':'Hedgehog',
'IL1':'Interleukin-1',
'Leptin':'Leptin',
'TWEAK':'TNF-related weak inducer of apoptosis',
'RAGE':'Advanced glycation end-products ',
'IL-11':'Interleukin-11',
'BDNF':'Brain-derived neurotrophic factor',
'IL-7':'Interleukin-7',
'Oncostatin_M':'Oncostatin M',
'Alpha6Beta4Integrin':'Alpha6 Beta4 Integrin',
'Notch':'Notch',
'IL2':'Interleukin-2',
'BCR':'B cell receptor',
'TGF_beta_Receptor':'Transforming growth factor beta receptor',
'IL9':'Interleukin-9',
'IL4':'Interleukin-4',
'TCR':'T Cell Receptor',
'EGFR1':'Epidermal growth factor receptor',
'TNFalpha':'Tumor necrosis factor alpha',
'RANKL':'Receptor activator of nuclear factor kappa-B ligand',
'CRH':'Corticotropin-releasing hormone',
'IL3':'Interleukin-3',
'TSH':'Thyroid-stimulating hormone',
'IL5':'Interleukin-5',
'Wnt':'Wnt',
}

ABBREVIATED_NAMES = {'TSLP':'TSLP',
'ID':'ID',
'KitReceptor':'Kit',
'IL6':'IL6',
'AndrogenReceptor':'AR',
'FSH':'FSH',
'Prolactin':'Prolactin',
'Hedgehog':'Hedgehog',
'IL1':'IL1',
'Leptin':'Leptin',
'TWEAK':'TWEAK',
'RAGE':'Advanced glycation end-products ',
'IL-11':'IL11',
'BDNF':'BDNF',
'IL-7':'IL7',
'Oncostatin_M':'OSM',
'Alpha6Beta4Integrin':'a6b4Integrin',
'Notch':'Notch',
'IL2':'IL2',
'BCR':'BCR',
'TGF_beta_Receptor':'TGFb',
'IL9':'IL9',
'IL4':'IL4',
'TCR':'TCR',
'EGFR1':'EGFR1',
'TNFalpha':'TNFa',
'RANKL':'RANKL',
'CRH':'CRH',
'IL3':'IL3',
'TSH':'TSH',
'IL5':'IL5',
'Wnt':'Wnt',
}

def make_nodes(df: pd.DataFrame) -> set:
    """
    :df      pandas dataframe
    :returns set of nodes
    """
    ## ignore nodes SRC or SINK
    n1 = {(x[0],x[2]) for x in df.values if x[0] != 'SRC'}
    n2 = {(x[1],x[2]) for x in df.values if x[1] != 'SNK'}
    return n1.union(n2)

def make_edges(df: pd.DataFrame,directed=True) -> set:
    """
    :df      pandas dataframe
    :returns set of edge tuples
    """
    ## ignore edges incident to SRC or SINK
    if directed:
    	return {(x[0],x[1]) for x in df.values if x[0] != 'SRC' and x[1] != 'SNK'}
    else:
    	return {frozenset(x) for x in df.values if x[0] != 'SRC' and x[1] != 'SNK'}

def load_df_tab(name:str):
    return pd.read_csv(name,sep='\t',engine='python')

def get_ground_truth(fname):
	ground = load_df_tab(fname)
	nodes = make_nodes(ground[['#tail','head','pathway_name']])
	dir_edges = make_edges(ground[['#tail','head','pathway_name']],directed=True)
	undir_edges = make_edges(ground[['#tail','head','pathway_name']],directed=False)
	return nodes,dir_edges,undir_edges

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

def main():
	files = glob.glob(PATHWAY_DIR+'*-nodes.txt')
	filenames = set()
	for f in files:
		n = f.split('/')[-1].split('-')
		n = '-'.join(n[:-1])
		filenames.add(n)
		#print("'"+n+"':'',")
	for n in sorted(list(filenames),key=lambda x:FULL_NAMES[x]):
		nodefile = PATHWAY_DIR+'/'+n+'-nodes.txt'
		edgefile = PATHWAY_DIR+'/'+n+'-edges.txt'
		nodes,dir_edges,undir_edges = get_ground_truth(edgefile)
		sources,sinks = get_labeled_nodes(nodefile)
		if len(sources) > 0 and len(sinks)>0:
			#print('%s & %s & %d & %d & %d & %d & %d\\\\' % (FULL_NAMES[n],ABBREVIATED_NAMES[n],len(nodes),len(sources),len(sinks),len(dir_edges),len(undir_edges)))
			print('%s & %d & %d & %d & %d & %d\\\\' % (ABBREVIATED_NAMES[n],len(nodes),len(sources),len(sinks),len(dir_edges),len(undir_edges)))

if __name__ == '__main__':
	main()