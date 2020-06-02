#
# Tobias Rubel | rubelato@reed.edu
#
# this program makes whisker plots of fmaxs for methods on pathways
import sys
import pandas as pd
import matplotlib.pyplot as plt
import os
DPATH = 'fmaxs'

def plot(pname: str) -> None:
    global DPATH
    pfile = os.path.join(DPATH,next(x for x in os.listdir(DPATH) if pname in x))
    dat = pd.read_csv(pfile,index_col=0)
    mthds = ['RWR_a0.85-t0.5','PCSF_r5-b1-w5-g3','SP','PL_k500','BTB','RN_y20']
    dat.columns = mthds
    dat.plot.box()
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('{}-boxplot.pdf'.format(pname))
    
def main(argv):
    pathway = argv[1]
    plot(pathway)

if __name__ == "__main__":
    main(sys.argv)
