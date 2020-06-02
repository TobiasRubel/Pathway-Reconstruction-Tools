## lists all pathways for human.

import collections
#from restful_lib import Connection
import requests ## udpate 2020
import sys

kegg_url = "http://rest.kegg.jp"
#conn = Connection(kegg_url)
org = sys.argv[1].lower()

REQUEST_ = '/list/pathway/' + org
org_pathways = requests.get(kegg_url+REQUEST_ , headers={'Accept':'text/json'})
org_dict = {}

myList = org_pathways.text.split('\n')
for l in myList:
    if '\t' in l:
        newl = l.split('\t')
        org_dict[newl[0]] = newl[1]

print '# %d %s pathways.'  % (len(org_dict.keys()) , org)
for p in org_dict.keys():
    print '%s --> %s' % (p,org_dict[p])

