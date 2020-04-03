#
# Tobias Rubel | rubelato@reed.edu
# Reed CompBio
# 
# just plots the data csvs. Probably not going to be used for much.
#
import sys
import pandas as pd
import matplotlib.pyplot as plt

def plot(df):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(list(df['labels']),list(df['edges']))
    plt.xlabel('number of sources and sinks')
    plt.ylabel('number of edge')
    plt.savefig('edge_vs_src+snk.pdf')


def main(argv):
    df = pd.read_csv(argv[1],index_col=0)
    plot(df)

if __name__=="__main__":
    main(sys.argv)



