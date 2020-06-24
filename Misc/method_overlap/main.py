#
# Tobias Rubel | rubelato@reed.edu
# Reed CompBio
#
# This program computes (and plots) the overlap between various methods.

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#global variables

DDIR = '../../Validation/PR/refactor-test-data'

#fetch data directories
COMPDIRS = [os.path.join(DDIR,x) for x in os.listdir(DDIR) if 'composite' in x and not ('Hybrid' in x or 'Perfect' in x or 'PRAUG' in x or 't0.5' in x)]
print(COMPDIRS)
pdict = {x.split('/')[-1].split('_')[0]:x for x in COMPDIRS}
print(pdict)

#load data into a dictionary
npdict = dict()
for p in pdict:
    df = pd.read_csv(os.path.join(pdict[p],'ranked-edges.csv'),sep='\t')
    npdict[p] = {x for x in df[['#tail','head']].stack()}

#define overlap function

sim = lambda x,y: len(x.intersection(y))/len(x)

#compute sims

simdict = {p:[sim(npdict[p],npdict[m]) for m in npdict] for p in npdict}

df = pd.DataFrame(simdict)
df.index = list(simdict.keys())
print(df)
sns.heatmap(df, annot=True,cmap=sns.cm.rocket_r)
locs,labels = plt.xticks()
sdict = {'PL':15403,'RWR':6074,'RN':791,'PCSF':677,'SP':2281,'BTB':788}
plt.xticks(locs,labels=['{}\n{}'.format(l,sdict[l]) for l in df.index])
plt.tight_layout()
plt.savefig('composite_method_overlap.pdf')
plt.show()
"""
plt.pcolor(df)
plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)
plt.show()
"""
