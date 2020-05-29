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

#DEST_PATH = '../Validation/PR/data'
DEST_PATH = '/Volumes/compbio/2020-05-PRAUG/runs'

#for each pathway, plots are to be deposited here.
#PLOT_PATH = '../Validation/PR/plots'
PLOT_PATH = '../plots'
#all input data is placed here...
DATA_PATH = '../Pathways'

#except the interactome...
INTERACTOMES = {
    '2018':'../Interactomes/PathLinker_2018_human-ppi-weighted-cap0_75.txt',
    '2015':'../Interactomes/background-interactome-pathlinker-2015.txt'
    }
INTERACTOME = None ## will be set in main.

METHOD_ABBREVIATIONS = {'run_PathLinker':'PL',
    'run_BowtieBuilder':'BTB',
    'run_ResponseNet':'RN',
    'run_PCSF':'PCSF',
    'run_RWR':'RWR',
    'run_ShortestPaths':'SP',
    }

EXAMPLE_CONFIG = '../Validation/PR/config.conf'
#after running an algorithm, we need to create a folder for computing PR

## TODO: global variables don't need to be specified in functions if they are accessed (only if they are WRITTEN/UPDATED)

def report(algorithm:str, interactome: str, pathway:str, args:argparse.Namespace) -> None:
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
    DEST = get_outdir(algorithm, interactome, pathway,args)
    DEST = os.path.abspath(DEST) ## to avoid relative path errors
    if not os.path.isdir(DEST):
        os.mkdir(DEST)
    print('report function: writing to DEST',DEST)

    #populate directory (make sure it's a csv file)
    #print([x for x in os.listdir('.') if prediction in x])
    pred = '{}.csv'.format(algorithm)
    edest = os.path.join(DEST,'ranked-edges.csv')
    shutil.move(pred,edest)

    #put pathway name in prediction file
    df = pd.read_csv(edest,sep='\t')
    df['pathway_name'] = [pathway for x in df.index]
    df.to_csv(edest,sep='\t',index=False)

    # relative links need to be three more levels deep (to be reference from Validation submodule)
    #exceptions occur when the symlink already exists
    groundfile = os.path.join(DATA_PATH,'{}-edges.txt'.format(pathway))
    groundfile = os.path.abspath(groundfile) ## to avoid relative path errors
    if not os.path.isfile(os.path.join(DEST,'ground.csv')):
        if not os.path.isfile(groundfile):
            sys.exit("ERROR: ground truth edges file doesn't exist:",groundfile)
        #print('SYMLINK',groundfile,os.path.join(DEST,'ground.csv'))
        os.symlink(groundfile,os.path.join(DEST,'ground.csv'))

    pathway_nodes = groundfile.replace('-edges','-nodes')
    pathway_nodes = os.path.abspath(pathway_nodes) ## to avoid relative path errors
    if not os.path.isfile(os.path.join(DEST,'ground-nodes.csv')):
        if not os.path.isfile(pathway_nodes):
            sys.exit("ERROR: ground truth nodes file doesn't exist:",pathway_nodes)
        #print('SYMLINK:',pathway_nodes,os.path.join(DEST,'ground-nodes.csv'))
        os.symlink(pathway_nodes,os.path.join(DEST,'ground-nodes.csv'))

    interactome = os.path.abspath(interactome) ## to avoid relative path errors
    if not os.path.isfile(os.path.join(DEST,'interactome.csv')):
        if not os.path.isfile(interactome):
            sys.exit("ERROR: interactome doesn't exist:",interactome)
        #print('SYMLINK:',interactome,os.path.join(DEST,'interactome.csv'))
        os.symlink(interactome,os.path.join(DEST,'interactome.csv'))

    shutil.copy(EXAMPLE_CONFIG,os.path.join(DEST,'config.conf'))

## return True if file exists, False otherwise.
def outfile_exists(algorithm:str, interactome: str, pathway:str,args:argparse.Namespace) -> str:
    DEST = get_outdir(algorithm,interactome,pathway,args)
    edge_dest = os.path.join(DEST,'ranked-edges.csv')
    return os.path.isfile(edge_dest)

def format_args(algorithm,args):
    if 'PL' in algorithm:
        return 'k{}'.format(args.k)
    elif 'RN' in algorithm:
        return 'y{}'.format(args.y)
    elif 'PCSF' in algorithm:
        return 'r{}-b{}-w{}-g{}'.format(args.r,args.b,args.w,args.g)
    elif 'RWR' in algorithm:
        return 'a{}-t{}'.format(args.a,args.t)
    else: # algorithm has no arguments. Return empty string
        return ''

# out directory format: ALG_INTERACTOME_PATHWAY_ARGS
def get_outdir(algorithm:str, interactome: str, pathway:str, args:argparse.Namespace) -> str:
    global DEST_PATH
    i = interactome.split('/')[-1].split('.')[0]
    suffix = format_args(algorithm,args)
    if suffix == '':
        name = '{}_{}_{}'.format(algorithm,args.interactome,pathway)
    else:
        name = '{}_{}_{}_{}'.format(algorithm,args.interactome,pathway,suffix)
    return os.path.join(DEST_PATH,name)

#
#standardize running methods via wrappers.
#if you want to add a method to this analysis, you should add
#a routine here.
#
#n.b. these are not expected to run on arbitrary machines without
#     some modification. Change the RUN assignments to the executable
#     file on your system.


## TODO: this makes an intermediate (.lp) file.  Where should these live?
def run_RWR(interactome:str, pathway:str, labeled_nodes:str, args:argparse.Namespace, subcall=False) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    method_name = METHOD_ABBREVIATIONS['run_RWR']
    if ( (not args.force) or (subcall and not args.force_subcall) ) and outfile_exists(method_name, interactome,pathway,args):
        print('{} Outfile Exists. Not overwriting: use --force to overwrite.'.format(method_name))
        return

    #set up what we need to execute
    RUN = 'Methods/RWR/rwr.py'
    verbose='True'
    CALL = 'python3 {} {} {} {} {} {}'.format(RUN,interactome,labeled_nodes,args.a,args.t,verbose)
    print(CALL)
    if not args.printonly:
        subprocess.call(CALL.split())
        report(method_name,interactome,pathway,args)

## TODO: this makes intermediate files.
## TODO: THIS IS AN UNDIRECTED GRAPH and thus the edges may be undirected.
def run_PCSF(interactome:str, pathway:str, labeled_nodes:str, args:argparse.Namespace, subcall=False) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    method_name = METHOD_ABBREVIATIONS['run_PCSF']
    if ( (not args.force) or (subcall and not args.force_subcall) ) and outfile_exists(method_name, interactome,pathway,args):
        print('{} Outfile Exists. Not overwriting: use --force to overwrite.'.format(method_name))
        return

    #set up what we need to execute
    RUN = 'Methods/PCSF/pcsf.py'
    verbose='True'
    CALL = 'python3 {} {} {} {} {} {} {} {}'.format(RUN,interactome,labeled_nodes,args.w,args.b,args.g,args.r,verbose)
    print(CALL)
    if not args.printonly:
        subprocess.call(CALL.split())
        report(method_name,interactome,pathway,args)

## TODO: this makes an intermediate (.lp) file.  Where should these live?
def run_ResponseNet(interactome:str, pathway:str, labeled_nodes:str, args:argparse.Namespace, subcall=False) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    method_name = METHOD_ABBREVIATIONS['run_ResponseNet']
    if ( (not args.force) or (subcall and not args.force_subcall) ) and outfile_exists(method_name, interactome,pathway,args):
        print('{} Outfile Exists. Not overwriting: use --force to overwrite.'.format(method_name))
        return

    #set up what we need to execute
    RUN = 'Methods/ResponseNet/response_net.py'
    verbose='True'
    CALL = 'python3 {} {} {} {} {}'.format(RUN,interactome,labeled_nodes,args.y,verbose)
    print(CALL)
    if not args.printonly:
        subprocess.call(CALL.split())
        report(method_name,interactome,pathway,args)

def run_BowtieBuilder(interactome:str, pathway:str, labeled_nodes:str, args:argparse.Namespace, subcall=False) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    method_name = METHOD_ABBREVIATIONS['run_BowtieBuilder']
    if ( (not args.force) or (subcall and not args.force_subcall) ) and outfile_exists(method_name, interactome,pathway,args):
        print('{} Outfile Exists. Not overwriting: use --force to overwrite.'.format(method_name))
        return

    #set up what we need to execute
    RUN = 'Methods/BowtieBuilder/bowtie_builder.py'
    verbose='True'
    CALL = 'python3 {} {} {} {}'.format(RUN,interactome,labeled_nodes,verbose)
    print(CALL)
    if not args.printonly:
        subprocess.call(CALL.split())
        report(method_name,interactome,pathway,args)

def run_ShortestPaths(interactome:str, pathway:str, labeled_nodes:str, args:argparse.Namespace, subcall=False) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """
    method_name = METHOD_ABBREVIATIONS['run_ShortestPaths']
    if ( (not args.force) or (subcall and not args.force_subcall) ) and outfile_exists(method_name, interactome,pathway,args):
        print('{} Outfile Exists. Not overwriting: use --force to overwrite.'.format(method_name))
        return

    #set up what we need to execute
    RUN = 'Methods/ShortestPaths/shortest_paths.py'
    verbose='True'
    CALL = 'python3 {} {} {} {}'.format(RUN,interactome,labeled_nodes,verbose)
    print(CALL)
    if not args.printonly:
        subprocess.call(CALL.split())
        report(method_name,interactome,pathway,args)

def run_PathLinker(interactome:str, pathway:str, labeled_nodes:str, args:argparse.Namespace, subcall=False) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """

    method_name = METHOD_ABBREVIATIONS['run_PathLinker']
    if ( (not args.force) or (subcall and not args.force_subcall) ) and outfile_exists(method_name, interactome,pathway,args):
        print('{} Outfile Exists. Not overwriting: use --force to overwrite.'.format(method_name))
        return

    #set up what we need to execute
    RUN = 'Methods/PathLinker/run.py'
    CALL = 'python2.7 {} -k {} -o PL {} {}'.format(RUN,args.k,interactome,labeled_nodes)
    print(CALL)
    if not args.printonly:
        subprocess.call(CALL.split()) # call function
        ## replace the outfile name with PL.csv
        shutil.move('PLk-{}-ranked-edges.txt'.format(args.k),'{}.csv'.format(method_name))
        report(method_name,interactome,pathway,args)

## TODO: Lost PerfectLinker (ground truth).  Add another option for algorithm here??
## TODO: right now DFS/weighted/edges are optional arguments - is there a better way to do this?
def PRAUG(algorithm:str, interactome:str, pathway:str, labeled_nodes:str, args:argparse.Namespace,DFS=True,weighted=False,edges=False) -> None:
    """
    :interactome   path/to/interactome
    :labeled_nodes path/to/source and dest node file
    :pathway       path/to/actual ground truth pathway
    :k             number of paths to compute (meaningless here)
    :returns       nothing
    :side-effect   makes a dest directory with predicted pathway
    """

    if DFS:
        method_name = 'PRAUG-{}'.format(algorithm)
    else:
        method_name = 'PRAUG-{}-BFS'.format(algorithm)
    if weighted:
        method_name +='-WEIGHTED'
    ## Note: assuming that if edges=True then the algorithm name reflects this (e.g. PRAUG-GT-EDGES)

    #print('DFS:',DFS,'WEIGHTED:',weighted,'EDGES',edges,'--> METHOD NAME',method_name)
    if not args.force and outfile_exists(method_name, interactome,pathway,args):
        DEST = get_outdir(method_name,interactome,pathway,args)
        edge_dest = os.path.join(DEST,'ranked-edges.csv')
        print('{} Outfile Exists: {}\n\tNot overwriting: use --force to overwrite.'.format(method_name,edge_dest))
        return

    if algorithm == 'PL':
        run_PathLinker(interactome,pathway,labeled_nodes,args,subcall=True)
    elif algorithm == 'BTB':
        run_BowtieBuilder(interactome,pathway,labeled_nodes,args,subcall=True)
    elif algorithm == 'RN':
        run_ResponseNet(interactome,pathway,labeled_nodes,args,subcall=True)
    elif algorithm == 'PCSF':
        run_PCSF(interactome,pathway,labeled_nodes,args,subcall=True)
    elif algorithm == 'RWR':
        run_RWR(interactome,pathway,labeled_nodes,args,subcall=True)
    elif algorithm == 'SP':
        run_ShortestPaths(interactome,pathway,labeled_nodes,args,subcall=True)
    elif algorithm == 'GT-NODES' or algorithm == 'GT-EDGES':
        predictions = labeled_nodes.replace('-nodes.txt','-edges.txt')
    else:
        sys.exit('ERROR: algorithm "%s" is not implemented.' % (algorithm))

    ## ground-truth nodes and edges already set the prediction file; for all others,
    ## get the ranked-edges.csv file.
    if algorithm != 'GT-NODES' and algorithm != 'GT-EDGES':
        predictions = os.path.join(get_outdir(algorithm,interactome,pathway,args),'ranked-edges.csv')

    #set up what we need to execute
    ## TODO: change name of PerfectLinker/PL.py to PRAUG

    if edges:
        run_type = 'edges'
    else:
        run_type = 'nodes'

    ## VARIANTS is a dictionary of dictionary: VARIANTS[DFS][weighted] willreturn a tuple of (run,outfile).
    VARIANTS = {True:{True:None,False:None},False:{True:None,False:None}}
    if algorithm != 'GT-NODES' and algorithm != 'GT-EDGES':
        VARIANTS[True][True] = ('Methods/PerfectLinker-DFS-Weighted/PL.py','{}-{}-PerfectLinker-DFS-Weighted.csv'.format(algorithm,run_type))
        VARIANTS[True][False] = ('Methods/PerfectLinker/PL.py','{}-{}-PerfectLinker.csv'.format(algorithm,run_type))
        VARIANTS[False][True] = ('Methods/PerfectLinker-BFS-Weighted/PL.py','{}-{}-PerfectLinker-BFS-Weighted.csv'.format(algorithm,run_type))
        VARIANTS[False][False] = ('Methods/PerfectLinker-BFS/PL.py','{}-{}-PerfectLinker-BFS.csv'.format(algorithm,run_type))
    else: # outfile changes for ground truth runs.
        VARIANTS[True][True] = ('Methods/PerfectLinker-DFS-Weighted/PL.py','Pathways-{}-PerfectLinker-DFS-Weighted.csv'.format(run_type))
        VARIANTS[True][False] = ('Methods/PerfectLinker/PL.py','Pathways-{}-PerfectLinker.csv'.format(run_type))
        VARIANTS[False][True] = ('Methods/PerfectLinker-BFS-Weighted/PL.py','Pathways-{}-PerfectLinker-BFS-Weighted.csv'.format(run_type))
        VARIANTS[False][False] = ('Methods/PerfectLinker-BFS/PL.py','Pathways-{}-PerfectLinker-BFS.csv'.format(run_type))
    RUN,outfile = VARIANTS[DFS][weighted]

    CALL = 'python3 {} {} {} {} {}'.format(RUN,run_type,interactome,predictions,labeled_nodes)
    print(CALL)
    if not args.printonly:
        subprocess.call(CALL.split())
        ## replace the outfile name with PL.csv (TODO this hsould happen in the method)
        shutil.move(outfile,'{}.csv'.format(method_name))
        report(method_name,interactome,pathway,args)
'''
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
'''

def fetch_arguments(pathways,args):
    global DATA_PATH, INTERACTOME
    arguments = []
    for pname in pathways:
        nodefile = os.path.join(DATA_PATH,'{}-nodes.txt'.format(pname))
        arguments.append([INTERACTOME,pname,nodefile,args])
    return arguments

def main(argv):
    global INTERACTOME, INTERACTOMES

    #initialize some values
    #with open('main.py','r') as f:
    #    ALL_METHODS = re.findall('run_.*(?=\()',f.read())[:-1]
    ALL_METHODS = ['run_PathLinker','run_ShortestPaths','run_RWR','run_ResponseNet','run_BowtieBuilder','run_PCSF']

    ## PATHWAYS directory includes allNPnodes.txt and NP_pathways.zip; make sure that '-' is in the variable.
    ALL_PATHWAYS = {x.split('-')[0] for x in os.listdir(DATA_PATH) if '-' in x and not 'all' in x}
    ALL_PATHWAYS.remove('ID')
    ALL_PATHWAYS.remove('RAGE')
    ALL_PATHWAYS.remove('IL11')

    ## parse arguments
    args,PATHWAYS,METHODS = parse_args(ALL_PATHWAYS,ALL_METHODS)

    ## set INTERACTOME
    INTERACTOME = INTERACTOMES[args.interactome]

    #####
    ## Run Pathway Reconstruction Methods
    if args.run or args.runpraug or args.upper_bounds:
        method_inputs = fetch_arguments(PATHWAYS, args)

        ## TODO: split up theses dependencies by writing directly to the output directory.
        ## only run the original methods if --run is specified but NOT --runpraug.
        ## If --runpraug is specified, then the original methods will be re-run as a sub-call.
        if args.run and not args.runpraug: ## run originl methods first.
            for arg in method_inputs:
                print('Running {} methods on pathway {}'.format(len(METHODS),arg[1]))
                lat = []
                for method in METHODS:
                    proc = Process(target=eval(method), args=(arg))
                    proc.start()
                    lat.append(proc)
                for l in lat:
                    l.join()

        if args.runpraug: ## run PRAUG methods next
            for arg in method_inputs:
                print('Running {} PRAUG-augmented methods on pathway {}'.format(len(METHODS),arg[1]))
                lat = []
                for method in METHODS:
                    abbrv = METHOD_ABBREVIATIONS[method]
                    proc = Process(target=eval('PRAUG'), args=([abbrv] + arg))
                    proc.start()
                    lat.append(proc)
                for l in lat:
                    l.join()

        if args.upper_bounds: ## run PRAUG-GT-Nodes and PRAUG-GT-Edges last.
            for arg in method_inputs:
                print('Running PRAUG Ground Truth methods on pathway {}'.format(arg[1]))
                lat = []
                proc = Process(target=eval('PRAUG'), args=(['GT-NODES'] + arg))
                proc.start()
                lat.append(proc)
                proc = Process(target=eval('PRAUG'), args=(['GT-EDGES'] + arg + [True,False,True]))
                proc.start()
                lat.append(proc)
                for l in lat:
                    l.join()

    #####
    ## Get all methods for computing precisioin and recall and/or plotting.
    if args.pr or args.plot:
        ORIG_AND_AUG_METHODS = []
        if args.run:
            ORIG_AND_AUG_METHODS += [METHOD_ABBREVIATIONS[m] for m in METHODS]
        if args.runpraug:
            ORIG_AND_AUG_METHODS += ['PRAUG-{}'.format(METHOD_ABBREVIATIONS[m]) for m in METHODS]
        if args.upper_bounds:
            ORIG_AND_AUG_METHODS += ['PRAUG-GT-NODES','PRAUG-GT-EDGES']

    #####
    ## Compute Precision/Recall
    if args.pr:
        RUN = '../Validation/PR/make_pr.py'
        print('Computing Precision and Recall for %d pathways and %d methods' % (len(PATHWAYS),len(METHODS)))
        for p in PATHWAYS:
            print('precision/recall: {}'.format(p))
            runs = " ".join([get_outdir(m,INTERACTOME,p,args).split('/')[-1] for m in ORIG_AND_AUG_METHODS])
            print('runs: {}'.format(runs))
            CALL = 'python3 {} {} {}'.format(RUN,DEST_PATH,runs)
            print(CALL)
            if not args.printonly:
                subprocess.call(CALL.split())

    #####
    ## Plot PR Curves
    if args.plot:
        RUN = '../Validation/PR/plot_pr.py'
        print('Plotting Precision and Recall for %d pathways and %d methods' % (len(PATHWAYS),len(METHODS)))
        for p in PATHWAYS:
            print('plotting precision/recall: {}'.format(p))
            runs = " ".join([get_outdir(m,INTERACTOME,p,args).split('/')[-1] for m in ORIG_AND_AUG_METHODS])
            CALL = 'python3 {} {} {} {}'.format(RUN,DEST_PATH,PLOT_PATH,runs)
            print(CALL)
            if not args.printonly:
                subprocess.call(CALL.split())

    #####
    ## Post GraphSpace Graphs
    if args.gs:
        print('WARNING: not refactored.\n')
        username,password = args.post_graphs
        IS_DRAFT=True
        for p in PATHWAYS:
            print(p)
            for m in METHODS:
                print('  ',m)
                algorithm = m.split('_')[-1]
                name = '%s-%s' % (p,algorithm)
                directory = get_outdir(algorithm,INTERACTOME,p,args.k)
                #print(username,password,name,directory,draft)
                RUN = '../Misc/visualize_networks/post_graphspace_graph.py'
                CALL = 'python3 {} {} {} {} {} {}'.format(RUN,username,password,name,directory,IS_DRAFT)
                print(CALL)
                subprocess.call(CALL.split())

    ### SPECIAL ACTIONS
    #####
    ## Compute Node Precision/Recall (Motivation figure)
    if args.node_pr:
        if METHODS != ALL_METHODS:
            print('WARNING: Running node motivation pr code for all methods, not -m methods.')
        all_methods = [METHOD_ABBREVIATIONS[m] for m in ALL_METHODS]
        RUN = '../Validation/PR/make_node_motivation_pr.py'
        for p in PATHWAYS:
            print('node motivation: {}'.format(p))
            runs = " ".join([get_outdir(m,INTERACTOME,p,args).split('/')[-1] for m in all_methods])
            CALL = 'python3 {} {} {} {}'.format(RUN,DEST_PATH,False,runs)
            print(CALL)
            if not args.printonly:
                subprocess.call(CALL.split())

            CALL = 'python3 {} {} {} {}'.format(RUN,DEST_PATH,True,runs)
            print(CALL)
            if not args.printonly:
                subprocess.call(CALL.split())
                print()

    #####
    ## Compute RWR parameter sweep fig
    if args.rwr_sweep:
        orig_tau = args.t
        possible_taus = [0.1,0.2,0.3,0.4,0.5,0.75]

        method_inputs = fetch_arguments(PATHWAYS, args)
        ## run PRAUG-RWR on each pathway
        for tau in possible_taus:
            args.t = tau
            for arg in method_inputs:
                print('RWR parameter sweep: {}'.format(arg[1]))
                ## running PRAUG will also run RWR if not already populated.
                PRAUG('RWR',*arg)

        ## compute precision/recall
        RUN = '../Validation/PR/make_pr.py'
        PLOT_RUN = '../Validation/PR/plot_pr.py'
        for p in PATHWAYS:
            print('RWR parameter sweep precision/recall: {}'.format(p))
            rwr_run = []
            for tau in possible_taus:
                args.t = tau
                rwr_run.append(get_outdir('PRAUG-RWR',INTERACTOME,p,args).split('/')[-1])
                if tau == possible_taus[-1]: # only plot the last RWR run (all previous ones are on top of each other).
                    rwr_run.append(get_outdir('RWR',INTERACTOME,p,args).split('/')[-1])
            runs = ' '.join(rwr_run)
            print('runs: {}'.format(runs))
            CALL = 'python3 {} {} {}'.format(RUN,DEST_PATH,runs)
            print(CALL)
            if not args.printonly:
                subprocess.call(CALL.split())
            CALL = 'python3 {} {} {} {}'.format(PLOT_RUN,DEST_PATH,PLOT_PATH,runs)
            print(CALL)
            if not args.printonly:
                subprocess.call(CALL.split())
        args.t = orig_tau ## reset tau

    #####
    ## Compute PL parameter sweep fig
    if args.pl_sweep:
        orig_k = args.k
        possible_ks = [50,100,500,1000,5000]

        method_inputs = fetch_arguments(PATHWAYS, args)
        ## run PRAUG-RWR on each pathway
        for k in possible_ks:
            args.k=k
            for arg in method_inputs:
                print('PL parameter sweep: {}'.format(arg[1]))
                ## running PRAUG will also run RWR if not already populated.
                PRAUG('PL',*arg)

        ## compute precision/recall
        RUN = '../Validation/PR/make_pr.py'
        PLOT_RUN = '../Validation/PR/plot_pr.py'
        for p in PATHWAYS:
            print('PL parameter sweep precision/recall: {}'.format(p))
            rwr_run = []
            for k in possible_ks:
                args.k=k
                rwr_run.append(get_outdir('PRAUG-PL',INTERACTOME,p,args).split('/')[-1])
                if k == possible_ks[-1]: # only plot the last PL run (all previous ones are on top of each other).
                    rwr_run.append(get_outdir('PL',INTERACTOME,p,args).split('/')[-1])
            runs = ' '.join(rwr_run)
            print('runs: {}'.format(runs))
            CALL = 'python3 {} {} {}'.format(RUN,DEST_PATH,runs)
            print(CALL)
            if not args.printonly:
                subprocess.call(CALL.split())
            CALL = 'python3 {} {} {} {}'.format(PLOT_RUN,DEST_PATH,PLOT_PATH,runs)
            print(CALL)
            if not args.printonly:
                subprocess.call(CALL.split())
        args.k = orig_k ## reset tau

    #####
    ## Compute benchmarks
    if args.benchmark:
        ## first run all benchmarks
        method_inputs = fetch_arguments(PATHWAYS,args)
        for arg in method_inputs:
            print('benchmark: {}'.format(arg[1]))
            lat = []
            for method in METHODS:
                abbrv = METHOD_ABBREVIATIONS[method]
                for DFS in [True,False]:
                    for WEIGHTED in [True,False]:
                        proc = Process(target=eval('PRAUG'), args=([abbrv] + arg + [DFS,WEIGHTED]))
                        proc.start()
                        lat.append(proc)
                        for l in lat:
                            l.join()

        ## compute precision/recall
        RUN = '../Validation/PR/make_pr.py'
        for m in METHODS:
            abbrvs = [METHOD_ABBREVIATIONS[m]]
            #benchmark_methods = [a for a in abbrvs]
            benchmark_methods = ['PRAUG-{}'.format(a) for a in abbrvs]
            benchmark_methods += ['PRAUG-{}-BFS'.format(a) for a in abbrvs]
            benchmark_methods += ['PRAUG-{}-WEIGHTED'.format(a) for a in abbrvs]
            benchmark_methods += ['PRAUG-{}-BFS-WEIGHTED'.format(a) for a in abbrvs]
            print('benchmark compute pr: %d pathways and %d methods' % (len(PATHWAYS),len(METHODS)))
            for p in PATHWAYS:
                print('making precision/recall: {}'.format(p))
                runs = " ".join([get_outdir(m,INTERACTOME,p,args).split('/')[-1] for m in benchmark_methods])
                print(' runs:',runs)
                CALL = 'python3 {} {} {}'.format(RUN,DEST_PATH,runs)
                print(CALL)
                if not args.printonly:
                    subprocess.call(CALL.split())

            ## plot precision/recall
            RUN = '../Validation/PR/plot_pr.py'
            print('benchmark plot pr: %d pathways and %d methods' % (len(PATHWAYS),len(METHODS)))
            for p in PATHWAYS:
                print('plotting precision/recall: {}'.format(p))
                runs = " ".join([get_outdir(m,INTERACTOME,p,args).split('/')[-1] for m in benchmark_methods])
                CALL = 'python3 {} {} {} {}'.format(RUN,DEST_PATH,PLOT_PATH,runs)
                print(CALL)
                if not args.printonly:
                    subprocess.call(CALL.split())


# argument parser
def parse_args(ALL_PATHWAYS,ALL_METHODS):
    parser = argparse.ArgumentParser()

    #add optional arguments
    parser.add_argument("--force",action="store_true",help="Run method even if outfile exists in the correct place. If --force and --runpraug are specified by --run is NOT specified, will force run PRAUG but won't overwrite original method runs.  Default False.")
    parser.add_argument("--printonly",action="store_true",help="Print the commands but do not run them. Default False.")
    parser.add_argument("-p","--pathways",metavar='STR',nargs="+",help="A list of pathways to make predictions for. Possible options are 'all' or:\n{}".format("\n".join(ALL_PATHWAYS)))
    parser.add_argument("-m","--methods",metavar='STR',nargs="+",help="A list of pathway reconstruction algorithms to run. Possible options are 'all' or:\n{}".format("\n".join(ALL_METHODS)))
    parser.add_argument('-i','--interactome',type=str,metavar='STR',default='2018',help='interactome version (either "2015" or "2018"). Default "2018".')

    group = parser.add_argument_group('Actions')
    group.add_argument("--run",action="store_true",help='Run original pathway reconstruction methods on -p and -m arguments.')
    group.add_argument("--runpraug",action="store_true",help='Run PRAUG-augmented pathway reconstruction methods on -p and -m arguments. PRAUG will also be included in PR and plotting.')
    group.add_argument("--upper_bounds",action="store_true",help='Include PRAUG_GT_Nodes and PRAUG_GT_Edges (GT: "ground truth"); formerly PerfectLinker-Nodes and PerfectLinker edges.')
    group.add_argument("--pr",action="store_true",help="Compute Precision/Recall plots for predictions based on -p and -m arguments.")
    group.add_argument("--plot",action="store_true",help="Plot Precision/Recall plots for predictions based on -p and -m arguments.")
    group.add_argument("--post_graphs",nargs=2,metavar='STR',default=[None,None],help="Post GraphSpace graphs for predictions based on -p and -m arguments. Takes TWO inputs: GraphSpace username and password.")

    group = parser.add_argument_group('Specialized Actions')
    group.add_argument("--node_pr",action="store_true",help='Plot node motivation precision/recall (which is a little different than the typical PR).')
    group.add_argument('--pl_sweep',action="store_true",help='Generate precision/recall for PL parameter sweep based on -p arguments (ignores any -m)')
    group.add_argument('--rwr_sweep',action="store_true",help='Generate precision/recall for RWR parameter sweep based on -p arguments (ignores any -m)')
    group.add_argument('--benchmark',action='store_true',help='Run DFS/BFS/weighted scenarios on methods and parameters.')

    group = parser.add_argument_group('Pathway Reconstruction Method Arguments')
    group.add_argument('-k',type=int,metavar='INT',default=500,help="PathLinker: number of shortest paths (k). Default 500.")
    group.add_argument('-y',type=int,metavar='INT',default=20,help="ResponseNet: sparsity parameter (gamma). Default 20.")
    group.add_argument('-r',type=int,metavar='INT',default=5,help="PCSF: terminal prize (rho). Default 5.")
    group.add_argument('-b',type=int,metavar='INT',default=1,help="PCSF: edge reliability (b). Default 1.")
    group.add_argument('-w',type=int,metavar='INT',default=5,help="PCSF: dummy edge weight (omega). Default 5.")
    group.add_argument('-g',type=int,metavar='INT',default=3,help="PCSF: degree penalty (g). Default 3.")
    group.add_argument('-a',type=float,metavar='FLOAT',default=0.85,help="RWR: teleportation probability (alpha). Default 0.85.")
    group.add_argument('-t',type=float,metavar='FLOAT',default=0.5,help="RWR: flux threshold (tau). Default 0.5.")

    args = parser.parse_args()

    ## check that at least one action is specified
    if not (args.run or args.runpraug or args.pr or args.plot or args.post_graphs[0] or args.node_pr or args.rwr_sweep or args.pl_sweep or args.benchmark):
        parser.print_help()
        sys.exit('\nERROR: At least one action must be specified. Exiting.')

    ## check that interactome is valid
    if args.interactome not in ['2015','2018']:
        parser.print_help()
        sys.exit('\nEROOR: --interactome must either be "2015" or "2018". Exiting.')

    ## check that at least one pathway is specified.
    ## (note this is OK if node_pr is True)
    if args.pathways == None and not args.node_pr:
        parser.print_help()
        sys.exit('\nERROR: At least one pathway (-p) must be specified. Exiting.')
    elif args.pathways == ['all']:
        PATHWAYS = ALL_PATHWAYS
    else:
        PATHWAYS = args.pathways
        if any([p not in ALL_PATHWAYS for p in PATHWAYS]):
            sys.exit('\nERROR: Pathways can be "all" or from {}'.format(ALL_PATHWAYS))

    ## check that at least one method is specified
    ## (note this is OK if node_pr is True)
    if args.methods == None and not (args.node_pr or args.rwr_sweep or args.pl_sweep):
        parser.print_help()
        sys.exit('\nERROR: At least one method (-m) must be specified. Exiting.')
    elif args.methods == ['all']:
        METHODS = ALL_METHODS
    elif args.methods == None:
        if not (args.node_pr or args.rwr_sweep or args.pl_sweep):
            sys.exit('\nERROR: Method must be specified for any action other than --node_pr and --rwr_sweep or --pl_sweep. Exiting.')
        METHODS = None
    else:
        METHODS = args.methods
        if any([m not in ALL_METHODS for m in METHODS]):
            sys.exit('\nERROR: Methods can be "all" or from {}'.format(ALL_METHODS))


    ## check that specified arguments are actually run. If not emit a warning.
    if args.k != 500 and not ('run_PathLinker' in METHODS or 'run_HybridLinker' in METHODS):
        print('WARNING: -k option is specified but PL or PRAUG-PL is not specified.')
    if args.y != 20 and not ('run_ResponseNet' in METHODS or 'run_HybridLinker_RN' in METHODS):
        print('WARNING: -y option is specified but RN or PRAUG-RN is not specified.')
    if (args.r != 5 or args.b != 1 or args.w != 5 or args.g != 3) and not ('run_PCSF' in METHODS or 'run_HybridLinker_PCSF' in METHODS):
        print('WARNING: at least one of -r, -b, -w, -g options are specified but PCSF or PRAUG-PCSF are not specified.')
    if (args.a != 0.85 or args.t != 0.5) and not ('run_RWR' in METHODS or 'run_HybridLinker_RWR' in METHODS):
        print('WARNING: at least one of -a and -t options are specified but RWR or PRAUG-RWR are not specified.')
    print()

    ## add a gs argument
    if args.post_graphs[0] != None:
        args.gs = True
    else:
        args.gs = False

    ## add force_subcall method.
    ## if --runpraug and --force are specified but --run is NOT specified, set force_subcall = False
    ## For example, when running PRAUG-PL, shoudl PL be --force run or not?
    if args.runpraug and args.force and not args.run:
        args.force_subcall = False
    elif args.force: ## --force specified, so set --force_subcall to be True
        args.force_subcall = True
    else: ## --force NOT specified, so set --force_subcall to be False
        args.force_subcall = False

    return args, PATHWAYS, METHODS

if __name__ == "__main__":
    main(sys.argv)
