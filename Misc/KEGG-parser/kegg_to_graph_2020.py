import glob
import sys

### UNIPROT MAPPING
import urllib.parse
import urllib.request

INTERACTOMES = {
    '2018':'../../Interactomes/PathLinker_2018_human-ppi-weighted-cap0_75.txt',
    '2015':'../../Interactomes/background-interactome-pathlinker-2015.txt'
    }
INTERACTOME = '2018'

PATHWAY_DIR = '../../KEGG-Pathways/'

## https://stackoverflow.com/questions/35569042/ssl-certificate-verify-failed-with-python3
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

def main():

	## get all nodes and undirected edges in interactome:
	nodes = set()
	edges = set()
	with open(INTERACTOMES[INTERACTOME]) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split()
			nodes.add(row[0])
			nodes.add(row[1])
			edges.add((row[0],row[1]))
			edges.add((row[1],row[0]))
	print('%d nodes and %d edges' % (len(nodes),len(edges)))
	
	files = glob.glob(PATHWAY_DIR+'*-entries.tsv')
	kegg_edges = {}
	for f in files:
		pathway = f.split('/')[-1].split('-')[0]
		if pathway == 'hsa04310':
			verbose=True
		else:
			verbose=False
			
		if verbose:
			print('FILE:',f)
		entries = {}
		t = 0
		to_map = set()
		with open(f) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				row = line.strip().split()
				entries[row[0]] = set()
				uids = set()
				for u in row[3].split('|'):
					uids.update(set(u.split(',')))
				
				for u in uids:
					if u in nodes:
						entries[row[0]].add(u)
						t+=1
					else:
						to_map.add(u)
		if verbose:
			print('  %d entries have %d orig uniprot ids from the interactome' % (len(entries),t))
		
		'''
		## map.
		if verbose:
			print('mapping entries that couldn\'t be found...')
		other_ids = known_by_another(to_map,nodes,verbose=verbose)

		t=0
		with open(f) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				row = line.strip().split()
				for u in row[3].split('|'):
					if u not in nodes and u in other_ids:
						for v in other_ids[u]:
							entries[row[0]].add(v)
							t+=1
		if verbose:
			print('  %d entries have %d mapped uniprot ids from the interactome' % (len(entries),t))
		'''
		f = f.replace('-entries','-relations')
		t=0
		ig = 0
		num_lines = 0
		missing_entries = 0
		with open(f) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				num_lines+=1
				row = line.strip().split()
				entryid1,entryid2,entrytype,interactiontype = row[:4]
				if not isIgnoredEdge(entryid1,entryid2,entrytype,interactiontype):
					if len(entries[entryid1])==0 or len(entries[entryid2])==0:
						missing_entries+=1
					if verbose:
						print('  adding entry %s-%s (%s) to (%s)' % ( entryid1,entryid2,' '.join(entries[entryid1]),' '.join(entries[entryid2])))
					## add it, no matter what direction it is.  
					for u in entries[entryid1]:
						for v in entries[entryid2]:
							if (u,v) in edges:
								if (u,v) not in kegg_edges:
									kegg_edges[(u,v)] = set()
									kegg_edges[(v,u)] = set()
								kegg_edges[(u,v)].add(pathway)
								kegg_edges[(v,u)].add(pathway)
								t+=1
								if verbose:
									print('  Edge (%s,%s) is ADDED' % (u,v))
							else:
								if verbose:
									print('  Edge (%s,%s) is SKIPPED (not in G).\n\t%s in G? %s\n\t%s in G? %s' % (u,v,u,u in nodes,v,v in nodes))
				else:
					ig+=1
					if verbose:
						print('  Edge (%s,%s) is IGNORED (type %s %s).' % (entryid1,entryid2,entrytype,interactiontype))
						for u in entries[entryid1]:
							for v in entries[entryid2]:
								print('\tcEdge (%s,%s) in G? %s\n\t%s in G? %s\n\t%s in G? %s' % (u,v,(u,v) in edges,u,u in nodes,v,v in nodes))
		if verbose:
			print('  %d edges are annotated as pathway %s; %d are ignored (planned) and are %d missing entries (unplanned). %d total lines read.' % (t,pathway,ig,missing_entries,num_lines))

	## write files
	out = open(PATHWAY_DIR+'/KEGG_Annotated_Interactions.txt','w')
	out.write('#u\tv\tpathways\n')
	for u,v in kegg_edges:
		out.write('%s\t%s\t%s\n' % (u,v, '|'.join([e for e in kegg_edges[(u,v)]])))
	out.close()
	print('wrote to '+PATHWAY_DIR+'/KEGG_Annotated_Interactions.txt')

### functions from original code:

######################################################################
## Determine whether edge is ignored from KEGG Edge Type and KEGG Edge SubType
## Edge information and quotes are from http://www.kegg.jp/kegg/xml/docs/
## Input: 
##    u,v (ints) internal KEGG Entry IDs.
##    entrytype (string) One of PPrel, ECrel, GErel, PCrel
##    interactiontypes (set of strings) a set of KEGG Subtypes (compound,activation,inhibition,etc)
## OUTPUT:
##    True if ignored, False otherwise
def isIgnoredEdge(u,v,entrytype,interactiontypes):
    
    ## DETERMINE WHETHER EDGE SHOULD BE IGNORED
    # GErel is gene expression interaction (TF -> Target Gene)
    # Ignore these edges
    if entrytype == 'GErel':
        return True

    # KEEP INDIRECT-EFFECT EDGES NOW
    # If the edge has an 'indirect-effect' subtype, ignore.
    # 'indirect-effect' = "indirect effect without molecular details"
    #elif 'indirect-effect' in interactiontypes:
    #    return True

    # If the edge has a 'state-change' subtype, ignore.
    elif 'state-change' in interactiontypes:
        #if u != v: # Verify that u=v
        #    print('WARNING! Edge with entry IDs %s - %s has "state-change" subtype' % (u,v))
        return True

    # If the edge has a 'missing-interaction' subtype, ignore.
    # 'missing-interaction' = missing interaction due to mutation, etc.
    elif 'missing-interaction' in interactiontypes:
        return True

    # If the edge has an 'expression' subtype, ignore.
    # FOR ALL HUMAN: 'expression' only appears alone or with 'compound'.
    # I hand-checked a few, and they all seem to be mislabeled as PPrel (should be GErel)
    elif 'expression' in interactiontypes:
        return True

    # A small handful (~10 for all human pathway entries) have NO subtype.  
    # Ignore these.
    elif interactiontypes=='None':
        return True

    # If we get this far, the edge is not ignored
    return False

def known_by_another(uids,nodes,verbose=False):
	if verbose:
		print('%d uids'%(len(uids)))

	url = 'https://www.uniprot.org/uploadlists/'
	
	## get common name for uid
	params = {
	'from': 'ACC+ID',
	'to': 'GENENAME',
	'format': 'tab',
	'query': ' '.join(list(uids))
	}

	data = urllib.parse.urlencode(params)
	data = data.encode('utf-8')
	req = urllib.request.Request(url, data)
	with urllib.request.urlopen(req) as f:
	   response = f.read()
	answer = response.decode('utf-8')

	common_names = {}
	common_names_rev = {}
	for line in answer.split('\n'):
		if 'From' in line: # skip header
			continue
		row = line.strip().split()
		if len(row)==2:
			if row[0] not in common_names:
				common_names[row[0]] = []
			if row[1] not in common_names_rev:
				common_names_rev[row[1]] = []
			common_names[row[0]].append(row[1])
			common_names_rev[row[1]].append(row[0])
	if verbose:
		print('%d have common names:'%(len(common_names)))

	all_common_names = set()
	for u in common_names.values():
		all_common_names.update(set(u))

	### get uids from name:
	## get common name for uid
	params = {
	'from': 'GENENAME',
	'to': 'ACC',
	'format': 'tab',
	'query': ' '.join(list(all_common_names))
	}

	data = urllib.parse.urlencode(params)
	data = data.encode('utf-8')
	req = urllib.request.Request(url, data)
	with urllib.request.urlopen(req) as f:
	   response = f.read()
	answer = response.decode('utf-8')

	mapped_ids = {}
	for line in answer.split('\n'):
		if 'From' in line: # skip header
			continue
		row = line.strip().split()
		if len(row) == 2:
			if row[1] in nodes:
				for kegg_uid in common_names_rev[row[0]]:				
					if kegg_uid not in mapped_ids:
						mapped_ids[kegg_uid]=set()
					mapped_ids[kegg_uid].add(row[1])
	if verbose:
		print('%d mapped ids' % len(mapped_ids))
	return mapped_ids

if __name__ == '__main__':
	main()

