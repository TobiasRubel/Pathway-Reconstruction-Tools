## ResponseNet implementation
## 
## From the paper: "Briding High-Throughput Genetic and Transcriptional 
## Data Reveals Cellular Reponses to Alpha-Synuclein Toxicity" 
## Yege-Lotem et. al., Nat. Genetics 2009
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2733244/
##
## This algorithm is used in the REsponseNet webserver - v3 last
## published Jul 2019:
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6602570/
##
## Modified from 2013 implementation by Anna Ritz

## This file contains the CPLEX version of the run() function.
## The instance of the problems are too large for CPLEX free license
## One can get an academic license, but I decided to reimplement this 
## using an open-source solver instead.

## Import Statements
import sys
import numpy.linalg
import random
from math import log as math_log
import time
import networkx as nx

#import cplex
#from cplex.exceptions import CplexError

## Static Variables
SUPERSOURCE='SUPERSOURCE'
SUPERSINK='SUPERSINK'
DELIM='_'

#########################################################
## run ResponseNet with a CPLEX solver. 
def run_cplex(G, sources, targets, gamma, outfile, force, verbose):

    flow_var_list = augment_graph(G, sources, targets,verbose)

    prob = cplex.Cplex()
    prob.objective.set_sense(prob.objective.sense.minimize)
    if verbose:
        print('\nInitializing datacheck parameter (only in verbose mode)')
        prob.parameters.read.datacheck.set(1)

    print('\n---- Setting Objective: %d variables' % (len(flow_var_list)))
    prob.variables.add(names=flow_var_list)
    if verbose:
        print('  Variable Names:', prob.variables.get_names()[:3],'...')

    ## for flow variables f_e:
    ##   min_f [ \sum_{e} -log(w_e) x f_e ] - [ \gamma x \sum_{SS-src edge} f_{SS-src}]
    ## = min_f [ \sum_{e \notin SS-src edges}  -log(w_e) x f_e ] + [ \sum_{SS-src edge} [ log(w_e) - \gamma ] x f_{SS-src} ]
    linear_tuples = []
    for u,v in G.edges():
        if u == SUPERSOURCE:
            linear_tuples.append((G[u][v]['flow_var'],G[u][v]['cost']-gamma))
        else:
            linear_tuples.append((G[u][v]['flow_var'],G[u][v]['cost']))
    prob.objective.set_linear(linear_tuples)
    if verbose:
        print('  Objective Coefficients', prob.objective.get_linear()[:3],'...')

    ## Counter for ALL CONSTRAINTS that will be written.
    constraintnum = 0
    
    ## First Constraint: Flow from Sources to Sinks must be balanced
    print('\n---- Setting Source/Sink Flow Constraints')
    ## variables are all SUPERSOURCE-source flow variables and all target-SUPERSINK flow variables.
    variables = [G[SUPERSOURCE][s]['flow_var'] for s in sources] + [G[t][SUPERSINK]['flow_var'] for t in targets]
    ## coefficients are +1 for SUPERSOURCE-source edges and -1 for target-SUPERSINK edges.
    coeff = [1.0]*len(sources) + [-1.0]*len(targets)
    eq = cplex.SparsePair(ind=variables,val=coeff)
    prob.linear_constraints.add(lin_expr=[eq],senses=['E'],rhs=[0],names=['c%d'%(constraintnum)])
    if verbose:
        print('  source/sink constraint','c%d' % (constraintnum),':',eq,'E 0')
    constraintnum += 1
    
    ## Next Constraints: Flow in = Flow out
    print('\n---- Setting Balance Constraints')
    count = 0
    for n in G.nodes():
        # skip source or sink node
        if n == SUPERSOURCE or n == SUPERSINK:
            continue
        ## variables are all SUCCESSORS of n and all PREDECESSORS of n.
        variables = [G[n][j]['flow_var'] for j in G.successors(n)] + [G[j][n]['flow_var'] for j in G.predecessors(n)]
        ## coefficients are +1 for SUCCESSORS of n and -1 for PREDECESSORS of n.
        coeff = [1.0]*G.out_degree(n) + [-1.0]*G.in_degree(n)
        eq = cplex.SparsePair(ind=variables,val=coeff)
        prob.linear_constraints.add(lin_expr=[eq],senses=['E'],rhs=[0],names=['c%d'%(constraintnum)])
        if verbose and count < 5:
            print('  node balance constraint','c%d' % (constraintnum),':',eq,'E 0')
        elif verbose and count == 5:
            print('...')
        count+=1
        constraintnum += 1

    ## Bounds: flow must be smaller than capacities
    print('\n---- Setting Bounds')
    count = 0
    for u,v in G.edges():
        prob.variables.set_lower_bounds(G[u][v]['flow_var'],0)
        prob.variables.set_upper_bounds(G[u][v]['flow_var'],G[u][v]['cap'])
        if verbose and count < 5:
            print('Bound: %.2f <= %s <= %.2f' % (0,G[u][v]['flow_var'],G[u][v]['cap']))
        elif verbose and count == 5:
            print('...')
        count+=1            

    ## Try solving.
    print('Solving',prob.problem_type[prob.get_problem_type()],'Problem...')
    outprefix = outfile.replace('.txt','')
    prob.write(outprefix+'_responsenet.lp')
    if verbose:
        print('Wrote LP file to ',outprefix+'_responsenet.lp')

    prob.solve()
    if verbose:
        print('Done Solving.')
        statID = prob.solution.get_status()
        print('Solution status = ' ,statID, ':',prob.solution.status[statID])
        print('Objective value  = ',prob.solution.get_objective_value())
        print('\nVariables:')
  
    ## get values
    for u,v in G.edges():
        G[u][v]['flow'] = prob.solution.get_values(G[u][v]['flow_var'])
        if verbose:
            print('%s: %.2f' % (G[u][v]['flow_var'],G[u][v]['flow']))

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

    ## remove all changes to the graph
    restore_graph(G,sources,targets,verbose)

    return

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


#########################################################
## Restores G in-place (no need to return this object)
def restore_graph(G,sources,targets,verbose):
    G.remove_node(SUPERSOURCE) # removes nodes and incident edges
    G.remove_node(SUPERSINK) # removes nodes and incident edges
    for u,v in G.edges():
        del G[u][v]['cap']
        del G[u][v]['flow_var']
        del G[u][v]['flow']
        del G[u][v]['index']
    return



