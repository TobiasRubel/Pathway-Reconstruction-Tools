#
# Tobias Rubel | rubelato@mongodb.reed.edu
# Reed CompBio
#
# this program plots a heatmap given a csv and a destination file.
# call it using python3 plot_heatmap.py [input] [output]
#
import sys
import pandas as pd
from matplotlib import pyplot
import seaborn as sns

d = None

def load_df(fname):
    return pd.read_csv(fname,index_col=0)

def plot(fname,sname):
    global d
    df = load_df(fname)
    d = df
    print(df)
    df = df.sort_values(by='HybridLinker',axis=1)
    print(df.loc['HybridLinker'])
    print(df)
    fig, ax = pyplot.subplots(figsize=(12,5))
    #plt = sns.heatmap(ax=ax,data=df, annot=False,square=True)
    plt = sns.clustermap(data=df, annot=False,square=False,metric="euclidean", method="ward")
    #pt = plt.get_figure()
    plt.savefig(sname)


def main(argv):
    fname,sname = argv[1],argv[2]
    plot(fname,sname)

if __name__ == "__main__":
    main(sys.argv)
