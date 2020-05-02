import networkx as nx
import sys
import os
import pandas as pd

## implements Bowtie Builder for directed, weighted graphs.
## Code modified from a version originally written by Allison Tegge @ VT
def run(G, sources, targets, verbose):
    all_sources_targets = sources.union(targets)

    ## pre-processing: get shortest paths from all sources.
    shortest_path_lens = {}
    for s in sources:
        if verbose:
            print('pre-processing',s)
        ## https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.shortest_paths.generic.shortest_path_length.html#networkx.algorithms.shortest_paths.generic.shortest_path_length
        shortest_path_lens[s] = nx.shortest_path_length(G,s,weight='cost')
        
    ## The quotes describing the steps are from the paper by Supper et al.

    ## "1. Initialize the pathway P with all nodes S \cap T..." (** should be \cup **)
    P = nx.DiGraph()
    P.add_nodes_from(sources) ## Add sourcse to P
    P.add_nodes_from(targets) ## Add targets to P

  
    ## "1 cont'd ...and flag all nodes in S \cap T as 'not visited'." (** should be \cup **)
    visited = set()
    not_visited = all_sources_targets

    ## "2. Calculate the distance matrix D_{|S|x|T|} between the
    ##  nodes in S and T with Dijkstra's algorithm."
    ## Note: this is negative log-transformed, so D is calculating the log of the objective function.
    best = get_best_pair(shortest_path_lens, sources, targets, visited, not_visited)
    
    ## "6. Repeat the steps 2-5 until every node in S is connected
    ##  to some node in T, and vice versa if such a path exists in G."
    i = 1
    while len(not_visited) > 0:
        if verbose:
            print(' Iteration %d: %d nodes visited and %d nodes not_visited' % (i,len(visited),len(not_visited)))

        ## "3. Select the shortest path in D that connects a 'not visited'
        ##  and a 'visited' node in P, or, if no such path exists, a 'not 
        ##  visited' node in S to a 'not visited' node in T."
        if best == None:
            if verbose:
                print(" Unvisited sources:", not_visited.intersection(sources))
                print(" Unvisted targets:", not_visited.intersection(targets))
            return P

        ## "4. Add the nodes and edges of the selected path to P..."
        path = nx.shortest_path(G,best[0],best[1],weight='cost')
        nx.add_path(P, path)
        for j in range(len(path)-1):
            if 'iter' not in P[path[j]][path[j+1]]: # if this is the first time we've seen this edge, note the iteration.
                P[path[j]][path[j+1]]['iter'] = i

        ##  "4 cont'd.  ...flag all nodes in the pathway as 'visited'."
        visited = visited.union(set(path))
        not_visited = all_sources_targets.difference(visited)
        
        ## "5. Update D to include all distances to the nodes in P^D(s,t)"
        new_sources = sources.union(P).difference(targets)
        new_targets = targets.union(P).difference(sources)
        
        ## update the shortest_path_lens dictionary with the new nodes.
        for s in new_sources.difference(shortest_path_lens.keys()):
            #if verbose:
            #    print(' pre-processing',s)
            shortest_path_lens[s] = nx.shortest_path_length(G,s,weight='cost')

        #if verbose:
        #    print(' getting best pair given ',new_sources,' sources and ',new_targets,' targets.')
        best = get_best_pair(shortest_path_lens, new_sources, new_targets, visited, not_visited)
        if verbose:
            print(' best:',best)

        ## "6. Repeat the steps 2–5 until every node in S is connected to some node in T, 
        ## and vice versa if such a path exists in G."
        i +=1

    return P


## gets the pair from sources to targets that have the smallest distance,
## preferring a (visited, not_visited) pair over a (not_visited,not_visited pair).
## This replaces Supper et al's "construct a matrix D representing all shortest paths".
##
## Requirements: shortests_paths contains all sources.
def get_best_pair(shortest_path_lens, sources, targets, visited, not_visited):
    best_one_visited = (None,None,1000000) ## best (s,t,dist) tuple if one is visited 
    best_neither_visited = (None,None,1000000) ## best (s,t,dist) tuple if neither is visited
    for s in sources:
        for t in targets:
            if s not in shortest_path_lens:
                sys.exit('ERROR: shortest_path_lens does not include source %s and/or target %s. Quitting.' % (s,t))
            if t not in shortest_path_lens[s]: # t is not eachable from s; skip
                continue
            path_len = shortest_path_lens[s][t]
            if (s in visited and t in not_visited) or (s in not_visited and t in visited) and path_len < best_one_visited[2]:
                best_one_visited = (s,t,path_len)
            elif s in not_visited and t in not_visited and path_len < best_neither_visited[2]:
                best_neither_visited = (s,t,path_len)

    ## prefer best_one_visited
    if best_one_visited[0] != None:
        best = [best_one_visited[0],best_one_visited[1]]
    ## otherwise, take best_neither_visited
    elif best_neither_visited[0] != None:
        best = [best_neither_visited[0],best_neither_visited[1]]
    else:
        print('Warning: no paths found that includes a not_visited node.')
        return None

    return best

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

def write_output(G,outfile,verbose=False):
    ## write the network
    out = open(outfile,'w')
    out.write('#tail\thead\titer_found\n')
    for u,v in G.edges():
        #print('%s\t%s\t%d' % (u,v,G[u][v]['iter']))
        out.write('%s\t%s\t%d\n' % (u,v,G[u][v]['iter']))
    out.close()
    if verbose:
        print('wrote to %s' % (outfile))
    return

def main(argv):
    interactome = argv[1]
    labeled_nodes = argv[2]
    verbose = bool(argv[3])

    print('preprocessing inputs...')
    interactome = df_to_graph(interactome)
    sources,sinks = get_labeled_nodes(labeled_nodes)
    ## sources/sinks need to be sets here.
    sources = set(sources)
    sinks = set(sinks)

    print('making prediction...')
    G = run(interactome, sources, sinks, verbose=verbose)

    print('saving prediction...')
    write_output(G,'bowtie_builder.csv',verbose=verbose)

    return

if __name__ == '__main__':
    main(sys.argv)
