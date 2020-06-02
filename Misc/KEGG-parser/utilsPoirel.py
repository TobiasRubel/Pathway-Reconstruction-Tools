'''
This is a modest collection of methods that I regularly reuse.
'''

import sys
import networkx as nx

def avg(list):
    sum=0
    for e in list:
        sum+=e
    return sum*1.0/len(list)

def readDict(f, fromCol=1, toCol=2, sep='\t'):
    '''
    Read the dict from the given tab-delimited file. The dict
    maps items in the fromCol to items in the toCol (1-based column index).
    '''
    itemMap = {}
    for line in open(f, 'r').readlines():
        if line=='':
            continue
        if line[0]=='#':
            continue
        items = line.rstrip().split(sep)
        if len(items)<max(fromCol, toCol):
            continue
        key = items[fromCol-1]
        val = items[toCol-1]
        if key=='':
            continue
        if val=='':
            continue
        itemMap[key] = val
    return itemMap

def readColumnsSep(f, sep='\t', *cols):
    '''
    Read multiple columns and return the items from those columns
    in each line as a tuple.

    foo.txt:
        a b c
        d e f
        g h i

    Calling "readColumnsSep('foo.txt', ' ',1, 3,)" will return:
        [(a, c), (d, f), (g, i)]

    '''
    if len(cols)==0:
        return []
    rows = []
    for line in open(f, 'r').readlines():
        if line=='':
            continue
        if line[0]=='#':
            continue
        items = line.rstrip().split(sep)
        if len(items)<max(cols):
            continue
        rows.append(tuple([items[c-1] for c in cols]))
    return rows

def readColumns(f, *cols):
    '''
    Read multiple columns and return the items from those columns
    in each line as a tuple.

    foo.txt:
        a b c
        d e f
        g h i

    Calling "readColumns('foo.txt', 1, 3)" will return:
        [(a, c), (d, f), (g, i)]

    '''
    if len(cols)==0:
        return []
    rows = []
    for line in open(f, 'r').readlines():
        if line=='':
            continue
        if line[0]=='#':
            continue
        items = line.rstrip().split('\t')
        if len(items)<max(cols):
            continue
        rows.append(tuple([items[c-1] for c in cols]))
    return rows



def readItemList(f, col=1, sep='\t'):
    '''
    Read the given column of the tab-delimited file f
    and return it as a list. Col is the 1-based column index.
    '''
    itemlist = []
    for line in open(f, 'r').readlines():
        if line=='':
            continue
        if line[0]=='#':
            continue
        items = line.rstrip().split(sep)
        if len(items)<col:
            continue
        itemlist.append(items[col-1])
    return itemlist

def readItemSet(f, col=1, sep='\t'):
    '''
    Read the given column of the tab-delimited file f
    and return it as a set. Col is the 1-based column index.

    A wrapper to readItemList, returning a set instead of a list.
    '''
    return set(readItemList(f, col=col, sep=sep))

def readEdges(netfile, sep='\t'):
    '''
    Read the set of edges in the given file. The file should have at least
    2 tab-delimited columns. A set of (col1, col2) pairs is read and returned.
    '''
    f = open(netfile, 'r')
    edges = set()
    for line in f:
        if line=='':
            continue
        if line[0]=='#':
            continue
        items = line.strip().split(sep)
        if len(items)<2:
            continue
        edges.add( (items[0], items[1]) )
    f.close()
    return edges


def readNetwork(edgeFile, name=''):
    '''
    Read an undirected network from the given tab-delimited file f.
    The file should have at lest two columns. Each line represents
    a single edge between the identifiers in columns 1 and 2. An optional
    3rd column may be provided with weights for each edge. All other
    columns are ignored.
    '''
    f = open(edgeFile, 'r')
    num_cols = len(f.readline().rstrip().split('\t'))
    f.close()
    g = nx.Graph()
    if num_cols<2:
        return None
    elif num_cols==2:
        g = nx.read_edgelist(edgeFile, comments='#', delimiter='\t')
        for u,v in g.edges():
            g[u][v]['weight'] = 1.0
    elif num_cols>2:
        edges = zip(readItemList(edgeFile,1), readItemList(edgeFile,2), [float(x) for x in readItemList(edgeFile,3)])
        for u,v,w in edges:
            g.add_edge(u, v, {'weight':w})
    g.name = name
    g.remove_edges_from(g.selfloop_edges())
    if '-' in g:
        g.remove_node('-')
    return g


def readNetworks(edgeFile):
    '''
    \CAUTION: YOU CAN ONLY USE THIS FUNCTION IF YOUR NETWORK-IDS ARE IN SORTED ORDER IN THE INPUT MYGRAPHS FILE
    \TODO: FIX THE ABOVE ISSUE
    Read a set of undirected networks from the given tab-delimited file f.
    The file should have at lest three columns. Each line represents
    a single edge in a graph (whose id is indicated by column 1)
    between the identifiers in columns 2 and 3.
    #(not yet implemented): An optional fourth column may be provided with weights for each edge. All other columns are ignored.
    '''
    f = open(edgeFile, 'r')
    num_cols = len(f.readline().rstrip().split('\t'))
    f.close()
    graphList = []  #list of graphs to return
    if num_cols != 3:
        return None
    elif num_cols==3:
        #g = nx.Graph()
        #g = nx.read_edgelist(edgeFile, comments='#', delimiter='\t')
        #for u,v in g.edges():
            #g[u][v]['weight'] = 1.0
    #elif num_cols>2:
        edges = zip(readItemList(edgeFile,1), readItemList(edgeFile,2), readItemList(edgeFile,3))
        prevName = ""   #we assume that there is no graph without name (i.e. 1st col is not empty)
        for name,u,v in edges:
            if name != prevName:
                g = nx.Graph()
                prevName = name
                g.name = name
                graphList.append(g)
            graphList[-1].add_edge(u, v, {'weight':1.00})

        for i in range(len(graphList)):
            graphList[i].remove_edges_from(graphList[i].selfloop_edges())
            if '-' in graphList[i]:
                graphList[i].remove_node('-')
    return graphList


def readDirectedNetwork(edgeFile, name=''):
    '''
    Read a directed network from the given tab-delimited file f.
    The file should have at lest two columns. Each line represents
    a single edge from the identifier in columns 1 to the id in column 2.
    An optional 3rd column may be provided with weights for each edge.
    All other columns are ignored.
    '''
    f = open(edgeFile, 'r')
    num_cols = len(f.readline().rstrip().split('\t'))
    f.close()
    g = nx.DiGraph()
    if num_cols<2:
        return None
    elif num_cols==2:
        g = nx.read_edgelist(edgeFile, comments='#', delimiter='\t',create_using=nx.DiGraph())
        for u,v in g.edges():
            g[u][v]['weight'] = 1.0
    elif num_cols>2:
        edges = zip(readItemList(edgeFile,1), readItemList(edgeFile,2), [float(x) for x in readItemList(edgeFile,3)])
        for u,v,w in edges:
            g.add_edge(u, v, {'weight':w})
    g.name = name
    g.remove_edges_from(g.selfloop_edges())
    if '-' in g:
        g.remove_node('-')
    return g


def jaccardIndex(A, B, step=1):
    '''
    Compute the Jaccard index between A and B at intervals of the
    given step size.

    Returns a list of tuples (x,y) where y is the JI at cutoff x.
    '''
    seenOnce = set()
    seenTwice = set()
    minRange = min(len(A), len(B))

    if step<1:
        step = 1

    #JIvalues stores a list of pairs [x,y], where 'y' is the JI at cutoff 'x'
    JIvalues = []
    for k in range(minRange):
        currA = A[k]
        currB = B[k]
        if currA in seenOnce:
            seenOnce.remove(currA)
            seenTwice.add(currA)
        elif currA not in seenTwice: #so currA hasn't been seen at all
            seenOnce.add(currA)
        if currB in seenOnce:
            seenOnce.remove(currB)
            seenTwice.add(currB)
        elif currB not in seenTwice: #so currB hasn't been seen at all
            seenOnce.add(currB)

        #get the various JI values
        if ((k+1)%step==0) or (k+1==minRange):
            currJI = 1.0*len(seenTwice)/(len(seenOnce)+len(seenTwice))
            JIvalues.append( (k+1, currJI) )
    return JIvalues


def generalKT(A, B, step=1):
    '''
    Compute the generalization fo Kentall's tau for top-k lists by Fagin, Kumar,
    and Sivakumar (http://epubs.siam.org/doi/abs/10.1137/S0895480102412856).
    Briefly, this method compute the number of mismatched pair ordering between
    the two lists.

    This method is explained nicely in Mccown and Nelson
    (http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.65.3888).
    '''
    mapA = {}
    mapB = {}
    for i,a in enumerate(A):
        mapA[a] = i
    for i,b in enumerate(B):
        mapB[b] = i

    minRange = min(len(A), len(B))

    if step<1:
        step = 1

    #vals stores a list of pairs [x,y], where 'y' is the generalized Kendall tau at cutoff 'x'
    vals = []
    mismatched = set()
    for k in range(minRange):
        a = A[k]
        b = B[k]

        rnkaInB = mapB.get(a, minRange+1)
        rnkbInA = mapA.get(b, minRange+1)

        for rnkTempbInB, tempb in enumerate(B[:k+1]):
            rnkTempbInA = mapA.get(tempb, minRange+1)
            if (k<rnkTempbInA and rnkaInB>rnkTempbInB) or (k>rnkTempbInA and rnkaInB<rnkTempbInB):
                mismatched.add( tuple(sorted((a, tempb))) )

        # only index to k here, since we compare a to b in the previous for loop
        # by indexing to k+1
        for rnkTempaInA, tempa in enumerate(A[:k]):
            rnkTempaInB = mapB.get(tempa, minRange+1)
            if (k<rnkTempaInB and rnkbInA>rnkTempaInA) or (k>rnkTempaInB and rnkbInA<rnkTempaInA):
                mismatched.add( tuple(sorted((tempa,b))) )

        #update vals
        if ((k+1)%step==0) or (k+1==minRange):
            numpairs = (k+1)*(k+1)
            currVal = 1.0-1.0*len(mismatched)/numpairs
            vals.append( (k+1, currVal) )
    return vals


def computePR(pos, neg, values,compressed=True):
    '''
    pos and neg are sets of positive and negative items, respectively, and
    values is a list of pairs (item, val) where items are ranked in increasing
    order of their val. Precision and recall are computed for every unique val
    in the list of (item,val) pairs, and PR values are returned as a list of
    pairs (precision, recall).

    New by Anna March 2015:  if compressed == False, then the precision and
    recall is recorded for every item in values (not unique).  This is useful
    for printing out outfiles of precision and recall.  Here, PR values are 
    returned as a list of (item,val,pos|neg|ignore,precision,recall) tuples.
    '''
    pr = []
    tp = 0.0
    fp = 0.0
    lastVal = None
    precision = 0
    recall = 0
    for item, val in sorted(values, key=lambda x: x[1]):
        if item in pos:
            tp += 1.0
        elif item in neg:
            fp += 1.0
        else:
            ## if compressed, then skip ignored.
            ## if not compressed, then add a line to PR.
            if not compressed:
                pr.append( (item,val,'ignore',precision,recall) )
            continue

        if tp==0 and fp==0:
            continue

        precision = tp/(tp+fp)
        recall = tp/len(pos)
        if compressed: 
            ## Chris's original code.
            if lastVal==None:
                pr.append( (precision, recall) )
            # if the val is the same as the lastVal, update 
            # the last (precision, recall) entry
            elif lastVal==val:
                pr[-1] = ( (precision, recall) )
            # only add a new entry if it is different from the previous
            elif pr[-1]!=(precision, recall):
                pr.append( (precision, recall) )
            lastVal = val
        else: 
            ## Anna's code. 
            ## write out (item,val,precision,recall).
            ## this DOES NOT "look back" in the list to replace ties.  
            ## We will do that last.
            if item in pos: # write item as a positive.
                pr.append( (item,val,'pos',precision,recall) )
            elif item in neg: # write item as a negative.
                pr.append( (item,val,'neg',precision,recall) )
            else:
                sys.exit('error: item has to be a pos or a neg by this point.')

    if not compressed:
        ## We have not yet addressed ties in the precision/recall list.
        ## when "compressed==True", then you only need to keep track of
        ## the last value (lastVal) and update accordingly.  Here, we need
        ## to go backwards through the list, and update row i with row i+1's
        ## precision and recall if row i's value is the same as row i+1's value.
        for i in reversed(range(len(pr)-1)):
            ## check if row i's value is the same as row i+1's
            if pr[i][1] == pr[i+1][1]:
                ## update with new tuple
                ## items 0-2 from pr[i], items 3-4 from pr[i+1]
                newrow = (pr[i][0],pr[i][1],pr[i][2],pr[i+1][3],pr[i+1][4])
                pr[i] = newrow
    return pr
