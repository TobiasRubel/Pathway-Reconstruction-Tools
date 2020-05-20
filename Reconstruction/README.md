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

Final file is here:

```
../../../plots/node-motivation-BowtieBuilder-PathLinker-ShortestPaths-PathLinker-2018-human-ppi-weighted-cap0-75-Wnt-10000.png
```

## 3. Benchmark figures:

Pick PL as input methods for all pathways; plot composite.

### 3a. DFS/BFS and weighted/unweighted (easy to do)

### 3b. varying k in PL (easy to do, take a bit of time run)

Lets set `k=[50,100,500,1000,5000,10000]`

### 3c. make variance plots (nearly done)

Pick pathways of interest (Wnt) & methods of interest.  

## 4. Generating Networks (Anna -- not started)

## Others?

### Prediction overlaps (venn or heatmap)
### subset composte on small/med pathways
