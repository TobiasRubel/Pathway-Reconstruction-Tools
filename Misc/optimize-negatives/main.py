#
# Tobias Rubel | rubelato@reed.edu
# Reed CompBio
#
# This program intelligently creates negative sets for each pathway, and also throws them into a composite negative set.
#
#
import os
import shutil
import sys
import pandas as pd
import subprocess
DATA_PATH = '../../Pathways'
PR_PATH = '../../Validation/PR/refactor-test-data'

def init():
    try:
        shutil.rmtree('negatives')
    except:
        pass
    os.mkdir('negatives')
def fmax(csvdoc: str) -> float:
   df = pd.read_csv(csvdoc)
   vs = [tuple(x) for x in df.values]
   f1 = lambda p,r:2*((p*r)/(p+r))
   return max([f1(*v) for v in vs])

def concat_negs(path: str) -> None:
    df = pd.DataFrame.from_csv(os.path.join(path,os.listdir(path)[0]))
    for x in os.listdir(path)[1:]:
        neg = pd.DataFrame.from_csv(os.path.join(path,x))
        df = pd.concat([df,neg])
    df.to_csv(os.path.join(path,'composite-negatives.csv'))

def main(argv):
    global DATA_PATH,PR_PATH
    n = int(argv[1])
    ALL_PATHWAYS = {x.split('-')[0] for x in os.listdir(DATA_PATH) if '-' in x and not 'all' in x}
    #ALL_PATHWAYS = {x for x in ALL_PATHWAYS if 'Wnt' in x}
    prgrm = 'python3 main.py'
    init()
    mthds = ['RWR_2018_{}_a0.85-t0.5','PCSF_2018_{}_r5-b1-w5-g3','SP_2018_{}','PL_2018_{}_k500','BTB_2018_{}','RN_2018_{}_y20']
    for P in ALL_PATHWAYS:
        try:
            mat = []
            for num in range(n):
                #compute PR
                cdir = os.getcwd()
                os.chdir('../../Reconstruction')
                CALL = '{} -p {} -m run_RWR run_PCSF run_PathLinker run_ResponseNet run_BowtieBuilder run_ShortestPaths --run --pr'.format(prgrm,P)
                print(CALL)
                subprocess.call(CALL.split())
                os.chdir(cdir)
                #temporarily save the negative set
                shutil.copy(os.path.join(PR_PATH,mthds[0].format(P),'negatives.csv'),'{}-{}-negatives.csv'.format(num,P))
                #compute fmax
                rn = [fmax(os.path.join(PR_PATH,m.format(P),'pr-edges.csv')) for m in mthds]
                mat.append(rn)
            #now we have a 6xn matrix in mat of fmax values. Let's treat it as a dataframe for ease
            df = pd.DataFrame(mat)
            df.to_csv('{}-fmax.csv'.format(P))
            print(df)
            #now lets modify it s.t. each value of the matrix is the distance from the median value
            #in the column. e.g.
            # [ 1 2 ]   [ 1 1 ]
            # [ 2 1 ] ->[ 0 0 ] 
            # [ 3 2 ]   [ 1 1 ]
            for col in df:
                med = df[col].median()
                df[col] = df[col].apply(lambda x: abs(x - med))
            print('transformed df\n\n{}'.format(df))
            #sum the rows of the dataframe to get a vector
            vec = [sum(df.loc[row]) for row in df.index]
            print('vec of medians:\n{}'.format(vec))
            #finally, take the min
            mn = 0
            for i in range(len(vec)):
                if vec[i] < vec[mn]:
                    mn = i
            #now we know that the ith negative set is the negative set we want. Let's save it.
            shutil.copy('{}-{}-negatives.csv'.format(mn,P),'negatives')
        except Exception as e:
            print('failed for pathway: {}'.format(P))
            print(e)

    concat_negs('negatives')
if __name__ == "__main__":
    main(sys.argv)




            
