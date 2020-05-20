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

DDIR = '../../Validation/PR/data'

#fetch data directories
COMPDIRS = [os.path.join(DDIR,x) for x in os.listdir(DDIR) if 'composit' in x and not ('Hybrid' in x or 'Perfect' in x or 'RWR' in x)]
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
sns.heatmap(df, annot=True)
plt.tight_layout()
plt.savefig('composite_method_overlap.pdf')
plt.show()
"""
plt.pcolor(df)
plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)
plt.show()
"""
