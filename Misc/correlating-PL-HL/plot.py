#
# Tobias Rubel | rubelato@reed.edu
# Reed CompBio
# 
# just plots the data csvs. Probably not going to be used for much.
#
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x,y = list(df['PL Nodes']),list(df['HL Edges'])
    ax.scatter(x,y)
    plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),color='r')
    plt.xlabel('$F_{max}$ of PathLinker Node Precision/Recall')
    plt.ylabel('$F_{max}$ of HybridLinker Edge Precision/Recall')
    plt.title('Let $K=500$, Pearson Correlation = {}'.format(df.corr()['PL Nodes']['HL Edges']))
    plt.savefig('correlation.pdf')


def main(argv):
    df = pd.read_csv(argv[1],index_col=0)
    plot(df)

if __name__=="__main__":
    main(sys.argv)



