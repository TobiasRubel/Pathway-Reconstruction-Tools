#
# Tobias Rubel | rubelato@reed.edu
# CompBio
#
# This script is designed to test various network reconstruction methods against
# one another on a large number of pathways.
#
# It is not designed to be especially flexible as a tool. To run it just
# execute python3 main.py.
#
# it makes various assumptions about directory structures, input formatting, etc.
#
# Depends: python3.5>=,...

import subprocess
import os
import sys
import pandas as pd
import shutil
import argparse
import re
import time
from multiprocessing import Process
#
#declare global variables
#

#for each algorithm/interactome/pathway, a directory named [algorithm]_[interactome]_[pathway]
#should be placed here
#DEST_PATH = '/home/tobias/Documents/Work/CompBio/PR/2018_interactome-data'

DEST_PATH = '../Validation/PR/data'
#for each pathway, plots are to be deposited here.
#PLOT_PATH = '../Validation/PR/plots'
PLOT_PATH = '../plots'
#all input data is placed here...
DATA_PATH = '../Pathways'

#except the interactome...
INTERACTOME = '../Interactomes/PathLinker_2018_human-ppi-weighted-cap0_75.txt'
#INTERACTOME = '../Interactomes/background-interactome-pathlinker-2015.txt'

EXAMPLE_CONFIG = '../Validation/PR/config.conf'
#after running an algorithm, we need to create a folder for computing PR

def report(algorithm:str, prediction: str,interactome: str, labeled_nodes: str,pathway:str, k: int) -> None:
    """
    :prediction
    :interactome
    :labeled_nodes
    :pathway
    :k
    :returns       nothing
    :side-effect   generates directory suitable for PR code
    """
    global DEST_PATH,EXAMPLE_CONFIG
    #make directory
    i = interactome.split('/')[-1].split('.')[0]
    p = '-'.join(labeled_nodes.split('/')[-1].split('-')[:-1])
    name = '_'.join([algorithm,i,p,str(k)])
    DEST = os.path.join(DEST_PATH,name)
    try:
        os.mkdir(DEST)
    except:
        print('directory "%s" already existed. Not overwriting.' % (DEST))
    #populate directory (make sure it's a csv file)
    #print([x for x in os.listdir('.') if prediction in x])
    pred = next(x for x in os.listdir('.') if prediction in x and '.csv' in x)
    edest = os.path.join(DEST,'ranked-edges.csv')
    os.replace(pred,edest)
    #put pathway name in prediction file
    df = pd.read_csv(edest,sep='\t')
    df['pathway_name'] = [pathway.split('/')[-1].split('-')[0] for x in df.index]
    df.to_csv(edest,sep='\t',index=False)

    #exceptions occur when the symlink already exists
    try:
        os.symlink(os.path.join('../../../',pathway),os.path.join(DEST,'ground.csv'))
    except Exception as e:
        pass
    try:
        pathway_nodes = pathway.replace('-edges','-nodes')
        os.symlink(os.path.join('../../../',pathway_nodes),os.path.join(DEST,'ground-nodes.csv'))
    except Exception as e:
        pass
    try:
        os.symlink(os.path.join('../../../',interactome),os.path.join(DEST,'interactome.csv'))
    except Exception as e:
        pass
    shutil.copy(EXAMPLE_CONFIG,os.path.join(DEST,'config.conf'))

## return True if file exists, False otherwise.
def outfile_exists(algorithm:str, prediction: str,interactome: str, labeled_nodes: str,pathway:str, k:int) -> str:
    DEST = get_outdir(algorithm,interactome,pathway,k)
    print(DEST)
    edge_dest = os.path.join(DEST,'ranked-edges.csv')
    return os.path.isfile(edge_dest)

def get_outdir(algorithm:str, interactome: str, pathway:str, k:int) -> str:
    global DEST_PATH
    i = interactome.split('/')[-1].split('.')[0]
    if '-' in pathway:
        p = '-'.join(pathway.split('/')[-1].split('-')[:-1])
    else:
        p = pathway
    name = '_'.join([algorithm,i,p,str(k)])
    DEST = os.path.join(DEST_PATH,name)
    return DEST

#
#standardize running methods via wrappers.
#if you want to add a method to this analysis, you should add
#a routine here.
#
#n.b. these are not expected to run on arbitrary machines without
#     some modification. Change the RUN assignments to the executable
#     file on your system.


## TODO: remove k as a requirement.
## TODO: HOW TO CHANGE ALPHA??? THis should be a parameter (like k).
## TODO: this makes an intermediate (.lp) file.  Where should these live?
def run_RWR(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    alpha = 0.85 ## FIX THIS FOR NOW -- this needs to change.
    thres = 0.50 ## default parameter now - predicted edes should capture this percentage of flux.
    if not force and outfile_exists('RWR', 'rwr_alpha{}_thres{}'.format(alpha,thres),interactome,labeled_nodes,pathway,k):
        print('ResponseNet Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/RWR/rwr.py'
    verbose='True'

    CALL = 'python3 {} {} {} {} {} {}'.format(RUN,interactome,labeled_nodes,alpha,thres,verbose)
    print(CALL)
    #execute script
    subprocess.call(CALL.split())
    report('RWR','rwr_alpha{}_thres{}'.format(alpha,thres),interactome,labeled_nodes,pathway,k)



## TODO: remove k as a requirement.
## TODO: HOW TO CHANGE GAMMA??? THis should be a parameter (like k).
## TODO: this makes intermediate files.
## TODO: THIS IS AN UNDIRECTED GRAPH and thus the edges may be undirected.
def run_PCSF(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    dummy_edge_weight = 5 ## FIX THIS FOR NOW -- this needs to change.
    edge_reliability=1
    degree_penalty=3
    prize=5
    outprefix = 'pcsf-w%d-b%d-g%d-p%d' % (int(dummy_edge_weight),int(edge_reliability),int(degree_penalty),int(prize))
    if not force and outfile_exists('PCSF', outprefix,interactome,labeled_nodes,pathway,k):
        print('PCSF Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/PCSF/pcsf.py'
    verbose='True'


    CALL = 'python3 {} {} {} {} {} {} {} {}'.format(RUN,interactome,labeled_nodes,dummy_edge_weight,edge_reliability,degree_penalty,prize,verbose)
    #execute script
    subprocess.call(CALL.split())

    ## this needs to change (as does the outprefix in pcsf.py)
    print('outprefix:',outprefix)
    report('PCSF',outprefix,interactome,labeled_nodes,pathway,k)

## TODO: remove k as a requirement.
## TODO: HOW TO CHANGE GAMMA??? THis should be a parameter (like k).
## TODO: this makes an intermediate (.lp) file.  Where should these live?
def run_ResponseNet(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    gamma = 20 ## FIX THIS FOR NOW -- this needs to change.
    if not force and outfile_exists('ResponseNet', 'response_net_gamma_{}'.format(gamma),interactome,labeled_nodes,pathway,k):
        print('ResponseNet Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/ResponseNet/response_net.py'
    verbose='True'

    CALL = 'python3 {} {} {} {} {}'.format(RUN,interactome,labeled_nodes,gamma,verbose)
    #execute script
    subprocess.call(CALL.split())
    os.remove(next(x for x in os.listdir('.') if '.lp' in x))
    report('ResponseNet','response_net_gamma_{}'.format(gamma),interactome,labeled_nodes,pathway,k)

## TODO: remove k as a requirement.
def run_BowtieBuilder(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    if not force and outfile_exists('BowtieBuilder', 'bowtie_builder',interactome,labeled_nodes,pathway,k):
        print('BowtieBuilder Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/BowtieBuilder/bowtie_builder.py'
    verbose='True'
    CALL = 'python3 {} {} {} {}'.format(RUN,interactome,labeled_nodes,verbose)
    #execute script
    subprocess.call(CALL.split())
    report('BowtieBuilder','bowtie_builder',interactome,labeled_nodes,pathway,k)

## TODO: remove k as a requirement.
def run_ShortestPaths(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """

    if not force and outfile_exists('ShortestPaths', 'shortest_paths',interactome,labeled_nodes,pathway,k):
        print('ShortestPaths Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/ShortestPaths/shortest_paths.py'
    verbose='True'
    CALL = 'python3 {} {} {} {}'.format(RUN,interactome,labeled_nodes,verbose)
    #execute script
    subprocess.call(CALL.split())
    report('ShortestPaths','shortest_paths',interactome,labeled_nodes,pathway,k)

def run_PathLinker(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    if not force and outfile_exists('PathLinker', 'PathLinker',interactome,labeled_nodes,pathway,k):
        print('PathLinker Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/PathLinker/run.py'
    CALL = 'python2.7 {} -k {} -o PathLinker {} {}'.format(RUN,k,interactome,labeled_nodes)
    #execute script
    subprocess.call(CALL.split())
    report('PathLinker','PathLinker',interactome,labeled_nodes,pathway,k)

def run_LocPL(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    if not force and outfile_exists('Localized-PathLinker', 'Localized-PathLinker',interactome,labeled_nodes,pathway,k):
        print('Localized-PathLinker Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/localized-pathlinker/Loc_PL_run.py'
    PROTEIN_LOC = 'Methods/localized-pathlinker/Data/Protein_Localization_Scores.txt'
    COMPPI = 'Methods/localized-pathlinker/Data/comppi--interactions--tax_hsapiens_loc_all.txt'
    CALL = 'python2.7 {} -k {} -o LocalizedPathLinker {} {} {} {}'.format(RUN,k,interactome,labeled_nodes,COMPPI,PROTEIN_LOC)
    #execute script
    subprocess.call(CALL.split())
    #make directory for run.
    report('Localized-PathLinker','LocalizedPathLinker',interactome,labeled_nodes,pathway,k)

def run_PerfectLinker_edges(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """

    if not force and outfile_exists('PerfectLinker-edges', 'edges-PerfectLinker',interactome,labeled_nodes,pathway,k):
        print('PerfectLinker-edges Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/PerfectLinker/PL.py'
    CALL = 'python3 {} edges {} {} {}'.format(RUN,interactome,pathway,labeled_nodes)
    print('*'*20)
    print(CALL)
    #execute script
    subprocess.call(CALL.split())
    report('PerfectLinker-edges','edges-PerfectLinker',interactome,labeled_nodes,pathway,k)


def run_PerfectLinker_nodes(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """

    if not force and outfile_exists('PerfectLinker-nodes', 'nodes-PerfectLinker',interactome,labeled_nodes,pathway,k):
        print('PerfectLinker-nodes Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/PerfectLinker/PL.py'
    CALL = 'python3 {} nodes {} {} {}'.format(RUN,interactome,pathway,labeled_nodes)
    #execute script
    subprocess.call(CALL.split())
    report('PerfectLinker-nodes','nodes-PerfectLinker',interactome,labeled_nodes,pathway,k)

def run_HybridLinker(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """

    if not force and outfile_exists('HybridLinker', 'HybridLinker',interactome,labeled_nodes,pathway,k):
        print('HybridLinker Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/HybridLinker/main.py'
    CALL = 'python3 {} 500 {} {} {}'.format(RUN,interactome,labeled_nodes,pathway)
    #execute script
    print(CALL)
    subprocess.call(CALL.split())
    report('HybridLinker','HybridLinker',interactome,labeled_nodes,pathway,k)

def run_HybridLinker_SP(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """

    if not force and outfile_exists('HybridLinker-SP', 'HybridLinker-SP',interactome,labeled_nodes,pathway,k):
        print('HybridLinker-SP Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/ShortestPaths/HL.py'
    CALL = 'python3 {} {} {}'.format(RUN,interactome,labeled_nodes,)
    #execute script
    subprocess.call(CALL.split())
    report('HybridLinker-SP','HybridLinker-SP',interactome,labeled_nodes,pathway,k)

def run_HybridLinker_BTB(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    if not force and outfile_exists('HybridLinker-BTB', 'HybridLinker-BTB',interactome,labeled_nodes,pathway,k):
        print('HybridLinker-BTB Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/BowtieBuilder/HL.py'
    CALL = 'python3 {} {} {}'.format(RUN,interactome,labeled_nodes,)
    #execute script
    subprocess.call(CALL.split())
    report('HybridLinker-BTB','HybridLinker-BTB',interactome,labeled_nodes,pathway,k)

def run_HybridLinker_RN(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """

    if not force and outfile_exists('HybridLinker-RN', 'HybridLinker-RN',interactome,labeled_nodes,pathway,k):
        print('HybridLinker-RN Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/ResponseNet/HL.py'
    gamma = 20
    CALL = 'python3 {} {} {} {}'.format(RUN,interactome,labeled_nodes,gamma)
    #execute script
    subprocess.call(CALL.split())
    report('HybridLinker-RN','HybridLinker-RN',interactome,labeled_nodes,pathway,k)

def run_HybridLinker_RWR(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """

    if not force and outfile_exists('HybridLinker-RN', 'HybridLinker-RN',interactome,labeled_nodes,pathway,k):
        print('HybridLinker-RN Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/RWR/HL.py'
    alpha = 0.85
    CALL = 'python3 {} {} {} {}'.format(RUN,interactome,labeled_nodes,alpha)
    #execute script
    subprocess.call(CALL.split())
    report('HybridLinker-RWR','HybridLinker-RWR',interactome,labeled_nodes,pathway,k)


def run_HybridLinker_PCSF(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """

    if not force and outfile_exists('HybridLinker-PCSF', 'HybridLinker-PCSF',interactome,labeled_nodes,pathway,k):
        print('HybridLinker-PCSF Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/PCSF/HL.py'
    verbose='True'
    dummy_edge_weight = 5 ## FIX THIS FOR NOW -- this needs to change.
    edge_reliability=1
    degree_penalty=3
    prize=5
    CALL = 'python3 {} {} {} {} {} {} {} {}'.format(RUN,interactome,labeled_nodes,dummy_edge_weight,edge_reliability,degree_penalty,prize,verbose)
    #execute script
    subprocess.call(CALL.split())
    report('HybridLinker-PCSF','HybridLinker-PCSF',interactome,labeled_nodes,pathway,k)


def run_HybridLinker_BFS(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """

    if not force and outfile_exists('HybridLinker-BFS', 'HybridLinker-BFS',interactome,labeled_nodes,pathway,k):
        print('HybridLinker-BFS Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/HybridLinker-BFS/main.py'
    CALL = 'python3 {} 500 {} {} {}'.format(RUN,interactome,labeled_nodes,pathway)
    #execute script
    subprocess.call(CALL.split())
    report('HybridLinker-BFS','HybridLinker-BFS',interactome,labeled_nodes,pathway,k)

def run_HybridLinker_BFS_Weighted(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """

    if not force and outfile_exists('HybridLinker-BFS-Weighted', 'HybridLinker-BFS-Weighted',interactome,labeled_nodes,pathway,k):
        print('HybridLinker-BFS-Weighted Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/HybridLinker-BFS-Weighted/main.py'
    CALL = 'python3 {} 500 {} {} {}'.format(RUN,interactome,labeled_nodes,pathway)
    #execute script
    subprocess.call(CALL.split())
    report('HybridLinker-BFS-Weighted','HybridLinker-BFS-Weighted',interactome,labeled_nodes,pathway,k)

def run_HybridLinker_DFS_Weighted(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    if not force and outfile_exists('HybridLinker-DFS-Weighted', 'HybridLinker-DFS-Weighted',interactome,labeled_nodes,pathway,k):
        print('HybridLinker-DFS-Weighted Outfile Exists. Not overwriting: use --force to overwrite.')
        return

    #set up what we need to execute
    RUN = 'Methods/HybridLinker-DFS-Weighted/main.py'
    CALL = 'python3 {} 500 {} {} {}'.format(RUN,interactome,labeled_nodes,pathway)
    #execute script
    subprocess.call(CALL.split())
    report('HybridLinker-DFS-Weighted','HybridLinker-DFS-Weighted',interactome,labeled_nodes,pathway,k)

def run_HybridLinker_paramsweep(interactome:str,labeled_nodes:str,pathway:str,k: int,force: bool) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    #set up what we need to execute
    RUN = 'Methods/HybridLinker/main.py'
    for i in [50,100,500,1000,5000,10000]:
        CALL = 'python3 {} {} {} {} {}'.format(RUN,i,interactome,labeled_nodes,pathway)
        #execute script
        print(CALL)
        subprocess.call(CALL.split())
        report('HybridLinker-{}'.format(i),'HybridLinker',interactome,labeled_nodes,pathway,k)

#running the experiment

def fetch_arguments(k,pathways,force):
    global DATA_PATH
    global INTERACTOME
    pathways = {frozenset((x,y)) for x in os.listdir(DATA_PATH) for y in os.listdir(DATA_PATH) if x.split('-')[:-1] == y.split('-')[:-1] and x != y and x.split('-')[0] in pathways}
    sorted_pathways = [sorted(tuple(x),key = lambda x: x.split('-')[-1]) for x in pathways]
    processed_pathways = [(os.path.join(DATA_PATH,x),os.path.join(DATA_PATH,y)) for (x,y) in sorted_pathways]
    arguments = [(INTERACTOME,y,x,k,force) for (x,y) in processed_pathways]
    return arguments

def pr_all():
    print('COMPUTING PR FOR METHODS')
    global DATA_PATH
    global PLOT_PATH
    global DEST_PATH
    #hdir = os.listdir('.')
    #ndir = '/home/tobias/Documents/Work/CompBio/PR'
    #os.chdir(ndir)
    pathway_names = set(['-'.join(x.split('-')[:-1]) for x in os.listdir(DATA_PATH)])
    pathway_names.remove('')
    print(pathway_names)
    RUN = '../Validation/PR/make_pr.py'
    for p in pathway_names:
        print('p: {}'.format(p))
        runs = " ".join([x for x in os.listdir(DEST_PATH) if p in x])
        print('runs: {}'.format(runs))
        CALL = 'python3 {} {} {}'.format(RUN,DEST_PATH,runs)
        print(CALL)
        #subprocess.call(CALL.split())
    #os.chdir(hdir)


def pr_node_motivation(pathway_names):
    print('COMPUTING PR FOR NODE MOTIVATION')
    global DATA_PATH
    global PLOT_PATH
    global DEST_PATH
    #hdir = os.listdir('.')
    #ndir = '/home/tobias/Documents/Work/CompBio/PR'
    #os.chdir(ndir)
    print(pathway_names)
    methods = ['ShortestPaths_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000',\
        'PathLinker_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000',\
        'BowtieBuilder_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000',\
        'PCSF_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000',\
        'ResponseNet_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000',
        'RWR_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000']
    RUN = '../Validation/PR/make_node_motivation_pr.py'
    for p in pathway_names:
        print('p: {}'.format(p))
        runs = " ".join(methods)
        CALL = 'python3 {} {} {} {}'.format(RUN,DEST_PATH,False,runs)
        print(CALL)
        subprocess.call(CALL.split())
        CALL = 'python3 {} {} {} {}'.format(RUN,DEST_PATH,True,runs)
        print(CALL)
        subprocess.call(CALL.split())
        print()

    #os.chdir(hdir)

def plot_all():
    global DATA_PATH
    global PLOT_PATH
    global DEST_PATH
    ## need to check if '-'
    pathway_names = set(['-'.join(x.split('-')[:-1]) for x in os.listdir(DATA_PATH) if '-' in x])
    print('pathwaynames:',pathway_names)
    RUN = '../Validation/PR/plot_pr.py'

    for p in pathway_names:
        runs = " ".join([x for x in os.listdir(DEST_PATH) if p in x])
        CALL = 'python3 {} {} {} {}'.format(RUN,DEST_PATH,PLOT_PATH,runs)
        print('CALL: {}'.format(CALL))
        subprocess.call(CALL.split())

def main(argv):
    #initialize some values
    with open('main.py','r') as f:
        METHODS = re.findall('run_.*(?=\()',f.read())[:-1]
    ## PATHWAYS directory includes allNPnodes.txt and NP_pathways.zip; make sure that '-' is in the variable.
    PATHWAYS = {x.split('-')[0] for x in os.listdir(DATA_PATH) if '-' in x and not 'all' in x}
    parser = argparse.ArgumentParser()
    #add optional arguments
    parser.add_argument("--force",action="store_true",help="Run method even if outfile exists in the correct place. Default false.")
    parser.add_argument("--pr_all",action="store_true",help="Compute Precision/Recall plots for all data in DEST_PATH")
    parser.add_argument("--plot_all",action="store_true",help="Plot Precision/Recall plots for all data in DEST_PATH")
    parser.add_argument("--node_pr",action="store_true",help='Plot node motiavation precision/recall, which is a little different than the typical PR.')
    parser.add_argument("-p","--pathways",nargs="+",help="A list of pathways to make predictions for. Possible options are:\n{}".format("\n".join(PATHWAYS.union({'all'}))))
    parser.add_argument("-m","--methods",nargs="+",help="A list of algorithms to run. Possible options are:\n{}".format("\n".join(METHODS+['all'])))
    parser.add_argument("--post_graphs",nargs=2,default=[None,None],help="Post GraphSpace graphs for predictions based on -p and -m arguments. Takes TWO inputs: GraphSpace username and password.")
    parser.add_argument('-k',type=int,default=10000,help="k value to pass to PathLinker. Defaults to 10,000.")
    parser.add_argument('-y',type=int,default=20,help="gamma value to pass to ResponseNet. Defaults to 20.")
    args = parser.parse_args()
    if args.pathways != ['all']:
        PATHWAYS = args.pathways
    if args.methods != ['all']:
        METHODS = args.methods
    else: ## remove paramsweep & locPL
        METHODS.remove('run_HybridLinker_paramsweep')
        METHODS.remove('run_LocPL')
    if METHODS != None and PATHWAYS != None:
        ARGS = fetch_arguments(args.k,PATHWAYS,args.force)
        for arg in ARGS:
                print('running all methods on {}'.format(arg))
                lat = []
                for method in METHODS:
                    proc = Process(target=eval(method), args=(arg))
                    proc.start()
                    lat.append(proc)
                for l in lat:
                    l.join()
    if args.pr_all:
        pr_all()
    if args.node_pr:
        pr_node_motivation(PATHWAYS)
    if args.plot_all:
        plot_all()
    if args.post_graphs[0] != None:
        # runs this serially.
        username,password = args.post_graphs
        is_draft=True
        for p in PATHWAYS:
            print(p)
            for m in METHODS:
                print('  ',m)
                algorithm = m.split('_')[-1]
                name = '%s-%s' % (p,algorithm)
                directory = get_outdir(algorithm,INTERACTOME,p,args.k)
                #print(username,password,name,directory,draft)
                RUN = '../Misc/visualize_networks/post_graphspace_graph.py'
                CALL = 'python3 {} {} {} {} {} {}'.format(RUN,username,password,name,directory,is_draft)
                print(CALL)
                subprocess.call(CALL.split())

if __name__ == "__main__":
    main(sys.argv)
