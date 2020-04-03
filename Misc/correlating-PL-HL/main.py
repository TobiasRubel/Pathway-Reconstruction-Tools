#
# Tobias Rubel | rubelato@reed.edu
# Reed CompBio
#
# this program correlates the success of PathLinker at predicting nodes with the success of
# HybridLinker at predicting edges.

#
# note to self (delete once finished):
# this needs to run PL and HL on the same pathways with the same k value (we can just fix it for now)
# and then compute node and edge PR. 
# then I need to record PL node PR vs HL edge PR in a csv
from helper import *
import re
import pandas as pd
print(DEST_PATH)
#for each algorithm/interactome/pathway, a directory named [algorithm]_[interactome]_[pathway]
#should be placed here
DEST_PATH = '/home/tobias/Documents/Work/CompBio/correlating-PL-HL/data'

print(DEST_PATH)
#for each pathway, plots are to be deposited here.
PLOT_PATH = '/home/tobias/Documents/Work/CompBio/PR/2018_interactome-plots'

#all input data is placed here...
DATA_PATH = '/home/tobias/Documents/Work/CompBio/ritz-data/pathways'

#except the interactome...
INTERACTOME = '/home/tobias/Documents/Work/CompBio/localized-pathlinker/Data/PathLinker_2018_human-ppi-weighted-cap0_75.txt'

#pathways
PATHWAYS = set([x.split('-')[0] for x in os.listdir('/home/tobias/Documents/Work/CompBio/ritz-data/pathways')])
PATHWAYS.remove('allNPnodes.txt')
PATHWAYS.remove('RAGE')
PATHWAYS.remove('ID')

#what k value to bake in
K = 500

#process a data directory to get the PR data

def fmax(csvdoc: str) -> float:
   df = pd.read_csv(csvdoc)
   vs = [tuple(x) for x in df.values]
   f1 = lambda p,r:2*((p*r)/(p+r))
   return max([f1(*v) for v in vs])

def associate_runs() -> list:
    global DEST_PATH,PATHWAYS
    fs = os.listdir(DEST_PATH)
    print(PATHWAYS)
    assoc = [tuple(sorted([x for x in fs if re.search("_{}_".format(p),x)])) for p in PATHWAYS]
    print(assoc)
    return assoc

def make_correlation() -> None:
    global DEST_PATH
    assoc = associate_runs()
    PL_Node = []
    HL_Edge = []
    for a in assoc:
        try:
            x,y = a
            PL_Node.append(fmax(os.path.join(DEST_PATH,os.path.join(y,'pr-nodes.csv'))))
            HL_Edge.append(fmax(os.path.join(DEST_PATH,os.path.join(x,'pr-edges.csv'))))
        except:
            print(a)
    df = pd.DataFrame({'PL Nodes':PL_Node,'HL Edges':HL_Edge})
    print(df)
    df.to_csv('correlation.csv')



        
    


def main():
    make_correlation()
    return
    global K
    #all the methods to use 
    METHODS = [run_PathLinker, run_HybridLinker]
    #METHODS = [run_HybridLinker_BFS_Weighted,run_HybridLinker_DFS_Weighted]
    ARGS = fetch_arguments(K)
    #run predictions for all pathways. 
    #n.b. we will try to start |METHODS| processes.
    for arg in ARGS:
        print('running all methods on {}'.format(arg))
        lat = []
        for method in METHODS:
            proc = Process(target=method, args=(arg))
            proc.start()
            lat.append(proc)
        for l in lat:
            l.join()
    pr_all()
    plot_all()
    


if __name__ == "__main__":
    main()
 

