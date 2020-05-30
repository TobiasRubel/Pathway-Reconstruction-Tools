## KEGG Pathway Parser

Parser originally written in 2013 at Virginia Tech (with TM Murali and Allison Tegge).  Re-run here using python2 after updating to use `requests` library. See the documentation in the doc/ directory.

### Get all pathways from KEGG

```
python get_interactions.py -s hsa -c uniprot -a -o ../../KEGG-Pathways/ > ../../KEGG-Pathways/hsa.out
```

### Get the list of pathways

```
python list_pathways.py hsa > ../../KEGG-Pathways/HSA_PATHWAY_LIST.txt
```

### Finally, convert the kegg entries and relations to edges.  

These are undirected and are only written if the edge exists in the interactome (currently hard-coded to the 2018 version).  This ignores edges like the original `kegg_to_graph.py` code, but simply outputs each undirected edge with a pipe-delimited list of pathway IDs that contain the interaction. 

```
python3 kegg_to_graph_2020.py 
```
