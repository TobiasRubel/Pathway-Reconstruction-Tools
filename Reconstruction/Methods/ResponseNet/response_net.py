## ResponseNet implementation
##
## From the paper: "Briding High-Throughput Genetic and Transcriptional
## Data Reveals Cellular Reponses to Alpha-Synuclein Toxicity"
## Yege-Lotem et. al., Nat. Genetics 2009
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2733244/
##
## This algorithm is used in the ResponseNet webserver - v3 last
## published Jul 2019:
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6602570/
##
## Modified from 2013 implementation by Anna Ritz

## Import Statements
import sys
import numpy.linalg
import random
from math import log as math_log
import networkx as nx
import pandas as pd


## WARNING. MIP doesn't work with Python <3.5.1 b/c type checking is not implemented.
from mip import *

## Static Variables
SUPERSOURCE='SUPERSOURCE'
SUPERSINK='SUPERSINK'
DELIM='_'

## TODO comment

#########################################################
## run ResponseNet with MIP.
## https://docs.python-mip.com/en/latest/
def run(G, sources, targets, gamma, lp_outfile, verbose=False):

    flow_var_list = augment_graph(G, sources, targets,verbose)

    m = Model(sense=MINIMIZE) # make a MIP model that will minimize

    print('\n---- Setting Objective: %d variables' % (len(flow_var_list)))
    weights_w = []
    values_x = []
    ## for flow variables f_e:
    ##   min_f [ \sum_{e} -log(w_e) x f_e ] - [ \gamma x \sum_{SS-src edge} f_{SS-src}]
    ## = min_f [ \sum_{e \notin SS-src edges}  -log(w_e) x f_e ] + [ \sum_{SS-src edge} [ log(w_e) - \gamma ] x f_{SS-src} ]
    for u,v in G.edges():
        G[u][v]['flow_var_ref'] = m.add_var(name=G[u][v]['flow_var'], lb=0, ub=G[u][v]['cap'])
        values_x.append(G[u][v]['flow_var_ref'])
        if u == SUPERSOURCE:
            weights_w.append(G[u][v]['cost']-gamma)
        else:
            weights_w.append(G[u][v]['cost'])
    m.objective = minimize(xsum(weights_w[i]*values_x[i] for i in range(len(values_x))))

    ## First Constraint: Flow from Sources to Sinks must be balanced
    print('\n---- Setting Source/Sink Flow Constraints')
    ## variables are all SUPERSOURCE-source flow variables and all target-SUPERSINK flow variables.
    variables = [G[SUPERSOURCE][s]['flow_var_ref'] for s in sources] + [G[t][SUPERSINK]['flow_var_ref'] for t in targets]
    ## coefficients are +1 for SUPERSOURCE-source edges and -1 for target-SUPERSINK edges.
    coeff = [1.0]*len(sources) + [-1.0]*len(targets)
    m += xsum(coeff[i]*variables[i] for i in range(len(variables))) == 0

    ## Next Constraints: Flow in = Flow out
    print('\n---- Setting Balance Constraints')
    for n in G.nodes():
        # skip source or sink node
        if n == SUPERSOURCE or n == SUPERSINK:
            continue
        ## variables are all SUCCESSORS of n and all PREDECESSORS of n.
        variables = [G[n][j]['flow_var_ref'] for j in G.successors(n)] + [G[j][n]['flow_var_ref'] for j in G.predecessors(n)]
        ## coefficients are +1 for SUCCESSORS of n and -1 for PREDECESSORS of n.
        coeff = [1.0]*G.out_degree(n) + [-1.0]*G.in_degree(n)
        m += xsum(coeff[i]*variables[i] for i in range(len(variables))) == 0

    ## write file
    m.write(lp_outfile)
    if verbose:
        print('Wrote LP file to ',lp_outfile)

    ## Try solving.
    ## See https://docs.python-mip.com/en/latest/quickstart.html#creating-models
    status = m.optimize()
    if verbose:
        if status == OptimizationStatus.OPTIMAL:
            print('optimal solution cost {} found'.format(m.objective_value))
        elif status == OptimizationStatus.FEASIBLE:
            print('sol.cost {} found, best possible: {}'.format(m.objective_value, m.objective_bound))
        elif status == OptimizationStatus.NO_SOLUTION_FOUND:
            print('no feasible solution found, lower bound is: {}'.format(m.objective_bound))

    ## set values & print if verbose is True
    for u,v in G.edges():
        G[u][v]['flow'] = G[u][v]['flow_var_ref'].x
        if G[u][v]['flow'] > 1e-6 and verbose:
            print('%s: %.2f' % (G[u][v]['flow_var'],G[u][v]['flow']))

    return G

#########################################################
## Augments G in-place (no need to return this object)
## Returns an ordered list of flow variables
def augment_graph(G,sources,targets,verbose):
    if verbose:
        print(' augmenting graph...')

    ## initialize edge variables with capacities, flow_var name, and flow.
    for u,v in G.edges():
        G[u][v]['cap'] = 1.0 ## edge capacity
        G[u][v]['flow_var'] = DELIM.join(['f',u,v]) ## flow variable name for ILP
        G[u][v]['flow'] = None

    ## (1) Add super source and super sink
    G.add_node(SUPERSOURCE)
    G.add_node(SUPERSINK)

    ## (2) Add dir. edges from supersource to receptors and tfs to supersink

    ## (2a) Add capacities to edges from supersource to receptors

    ## NOTE: capacities are uniform for both of these types of edges.
    ## The original paper changed this to the strength of
    ## each genetic hit.  PathLinker used uniform capcities.

    for s in sources:
        G.add_edge(SUPERSOURCE,s,weight=1/len(sources),cap=1/len(sources),flow_var=DELIM.join(['f',SUPERSOURCE,s]),flow=None)
        G[SUPERSOURCE][s]['cost'] = -math_log(max([0.000000001, G[SUPERSOURCE][s]['weight']]))/math_log(10)

    ## (2b) Add capacities to eges from tfs to supersink
    for t in targets:
        G.add_edge(t,SUPERSINK,weight=1/len(targets),cap=1/len(targets),flow_var=DELIM.join(['f',t,SUPERSINK]),flow=None)
        G[t][SUPERSINK]['cost'] = -math_log(max([0.000000001, G[t][SUPERSINK]['weight']]))/math_log(10)

    orderedvariables = [] # sets the order of flow variables
    i = 0
    for u,v in G.edges():
        G[u][v]['index'] = i
        orderedvariables.append(G[u][v]['flow_var'])
        i+=1

    if verbose:
        print('augmented graph has %d nodes and %d edges' % (G.number_of_nodes(), G.number_of_edges()))

    return orderedvariables

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
    edgefile = open(outfile, 'w')
    edgefile.write('#tail\thead\tflow\n')
    for u,v in G.edges():
        if u == SUPERSOURCE or v == SUPERSINK: # ignore edges from SUPERSOURCE or to SUPERSINK
            continue
        ## Do not write if G[u][v]['flow'] is 0
        if G[u][v]['flow'] != 0:
            edgefile.write('%s\t%s\t%.4f\n' % (u,v,G[u][v]['flow']))
    edgefile.close()
    if verbose:
        print('Wrote to ',outfile)
    return

def main(argv):
    interactome = argv[1]
    labeled_nodes = argv[2]
    gamma = float(argv[3])
    verbose = bool(argv[4])

    print('preprocessing inputs...')
    interactome = df_to_graph(interactome)
    sources,sinks = get_labeled_nodes(labeled_nodes)

    ## sources/sinks need to be sets here.
    sources = set(sources)
    sinks = set(sinks)

    print('making prediction...')
    G = run(interactome, sources, sinks, gamma, 'RN.lp', verbose=verbose)

    print('saving prediction...')
    write_output(G, 'RN.csv',verbose=verbose)

    return

if __name__ == '__main__':
    main(sys.argv)
