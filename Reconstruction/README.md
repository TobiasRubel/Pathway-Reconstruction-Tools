# Pipelines for Generating Figures

## 1. Generating the composite figures (nearly done except for RWR)

1. Run all PR and visualizations for existing methods

```
python3 main.py --pr_all --plot_all
````

2. Generate the pathway union

```
python3 gen_pathway_union.py
```

3. Make PR with composite files only

```
python3 make_pr.py data $(ls data | grep composit)
```

4. Plot PR with composite files only. First set `COMPOSITE=True` in `plot_pr.py`. Then, run

```
python3 plot_pr.py data plots $(ls data | grep composit)
```

## 2. Generating the node motivation figure (Anna -- a little broken)

1. Run the following code to get predictions for the five relevant methods.

```
python3 main.py -p Wnt -m run_PCSF run_ResponseNet run_BowtieBuilder run_ShortestPaths run_PathLinker run_RWR
```

2. Run the node PR for these, subsampling 50X nodes and also ignoring adjacents: (TODO need to add RN, PCSF)

```
python3 main.py --node_pr -p Wnt
```

3. Plot the methods. Set `COMPOSITE=False` and `NODE_MOTIVATION=True` in `plot_pr.py` and then run:

```
python3 ../Validation/PR/plot_pr.py ../Validation/PR/data ../../../plots PathLinker_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000 ShortestPaths_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000 BowtieBuilder_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000 RWR_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000 PCSF_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000 ResponseNet_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000
```

Final files are here:

```
../plots/node-motivation-BowtieBuilder-PCSF-PathLinker-RWR-ResponseNet-ShortestPaths-PathLinker-2018-human-ppi-weighted-cap0-75-Wnt-10000.pdng
../plots/node-motivation-BowtieBuilder-PCSF-PathLinker-RWR-ResponseNet-ShortestPaths-PathLinker-2018-human-ppi-weighted-cap0-75-Wnt-10000.pdf
```

## 3. Benchmark figures:

Pick PL as input methods for all pathways; plot composite.

### 3a. DFS/BFS and weighted/unweighted and PerfectLinker Nodes/Edges (easy to do)

### 3b. varying k in PL & PRAUG(PL) (easy to do, take a bit of time run)

Lets set `k=[50,100,500,1000,5000,10000]`

### 3c. make variance plots (nearly done)

Pick pathways of interest (Wnt) & methods of interest.  

## 4. Generating Networks

For any number of pathways and methods, we can generate two networks: (a) the subgraph of predictions, and (b) an induced subgraph of predictions on the ground truth nodes. To post graphs for `HybridLinker` on the `Wnt` pathway, run:

```
python3 main.py --post_graphs <USERNAME> <PASSWORD> -p Wnt -m run_HybridLinker
```

where `<USERNAME>` is your [GraphSpace](http://graphspace.org/) username and `<PASSWORD>` is your [GraphSpace](http://graphspace.org/) password.    To run all methods on the Wnt pathway, run: (TOD add BTB-HL)

```
python3 main.py --post_graphs <USERNAME> <PASSWORD> -p Wnt -m run_RWR run_PCSF run_ResponseNet run_BowtieBuilder run_ShortestPaths run_PathLinker run_HybridLinker run_HybridLinker_SP run_HybridLinker_RN run_HybridLinker_RWR run_HybridLinker_PCSF run_HybridLinker_BTB
```

When `IS_DRAFT=False` in `main.py`, the graphs are shared with the []'reconstruction-traversals' GraphSpace group](http://graphspace.org/groups/1268) (membership required).

## Others?

### Prediction overlaps (venn or heatmap)
### subset composte on small/med pathways

# Example Runs (for debuging/refactoring)

```
time nohup python3 main.py --run --runpraug --pr --plot -p all -m all 
 ```

 On my mac, saving to external SSD, the above command (with no pre-processed experiments) filled up 681M and took  

 ```
real	90m0.324s
user	220m55.093s
sys	7m33.662s
```
