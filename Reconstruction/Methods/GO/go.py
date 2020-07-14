#
# Tobias Rubel | rubelato@reed.edu
# Reed CompBio
# 
# This program imagines a list of go terms as a pathway reconstruction method.
# Given G,s,t the method returns a set of vertices, as well as the empty set for edges.
#
# This allows us to incorporate GO terms into our current formalism of the pathway
# reconstruction problem, but it is fairly impragmatic, since we don't want to have 
# to proliferate new methods for each pathway... Thus we may want to refactor it 
# to automatically figure out which list of GO terms to use based on the sources and 
# targets. As it stands I am using just Wnt terms.
#
import sys
import pandas as pd

def load_go_terms(csvname: str) -> pd.DataFrame:
    """
    :csvname
    :returns dataframe of Uniprot ids, gene name
    """
    df = pd.read_csv(csvname,names=['bioentity'])
    df['bioentity'] = df['bioentity'].apply(lambda x: x.split(':')[-1])
    return df

def write_go_terms(df: pd.DataFrame) -> None:
    ndf = pd.DataFrame({'#tail':df['bioentity'],'head':['NaN' for x in df['bioentity']]})
    ndf.to_csv('GO.csv',sep='\t')


def main(argv):
    write_go_terms(load_go_terms('Methods/GO/Wnt-go.txt'))


if __name__ == "__main__":
    main(sys.argv)
