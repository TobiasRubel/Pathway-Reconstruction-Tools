#
# Tobias Rubel | rubelato@mongodb.reed.edu
# Reed Compbio 
# 
# Given a set M of methods and a set P of pathways,
# this tool generates an |M|x|P| matrix of definite
# integrals of the precision-recall curves for 
# those methods on those pathways.
import os
import pandas as pd
import numpy as np
import re

from sklearn import metrics
from scipy import integrate

def load_pr(path: str,fname: str) -> pd.DataFrame:
    """
    :path    to csv
    :fname   of csv
    :returns DataFrame of precision,recall
    """
    return pd.read_csv(os.path.join(path,fname))

def AUC(df: pd.DataFrame) -> int:
    """
    :df      DataFrame of precision,recall
    :returns definite integral of df at max recall 
    """
    #sort the dataframe
    df = df.sort_values(by=['recall','precision'],ascending=[True,False])
    return metrics.auc(df['recall'],df['precision'])
    #return integrate.simps(df['precision'],df['recall'])

def process_method_pathway(method:str,pathway:str,path:str,k=10000) -> int:
    """
    :method  name of method
    :pathway name of pathway
    :path    to precision/recall file 
    :returns AUC(method on pathway)
    """
    r = '{}.*{}.*{}'.format(method,pathway,k)
    try:
        fname = next(y for y in [re.match(r,x) for x in os.listdir(path)] if y != None)[0]
        pth = os.path.join(path,fname)
        df = load_pr(pth,'pr.csv')
        return AUC(df)
    except Exception as e:
        print('failed to find file matching regex {}'.format(r))
        print(e)

def process_method(method: str,pathways: list, path: str,k=10000) -> dict:
    """
    :method   name of method
    :pathways names of pathways
    :path     to precision/recall file 
    :returns  {pathway: AUC(method on pathway) for pathway in pathways}
    """
    return {p: process_method_pathway(method,p,path,k) for p in pathways}

def generate_heatmap(methods: list,pathways: list, path: str, k = 10000) -> pd.DataFrame:
    """
    :methods  names of methods
    :pathways names of pathways
    :path     to precision/recall file 
    :returns  {pathway: AUC(method on pathway) for pathway in pathways}
    """
    df = pd.DataFrame()
    for method in methods:
        aucs = process_method(method,pathways,path,k)
        df[method] = aucs.values()
    df.index = aucs.keys()
    return df.T

#declare globals and define main

METHODS = ['HybridLinker',
           'PerfectLinker-nodes',
           'PerfectLinker-edges',
           'PathLinker',]

PATHWAYS = set([x.split('-')[0] for x in os.listdir('/home/tobias/Documents/Work/CompBio/ritz-data/pathways')])
PATHWAYS.remove('allNPnodes.txt')
PATHWAYS.remove('RAGE')
PATHWAYS.remove('ID')
print(PATHWAYS)

#Specify this.
DATA_PATH = '/home/tobias/Documents/Work/CompBio/PR/debug-data'

def main():
    global METHODS
    global PATHWAYS
    global DATA_PATH
    df = generate_heatmap(METHODS,PATHWAYS,DATA_PATH)
    df.to_csv('debug-test.csv')

if __name__ == "__main__":
    main()





    
