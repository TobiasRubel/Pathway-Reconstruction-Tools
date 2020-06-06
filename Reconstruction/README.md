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

## 2. Generating the node motivation figure

Run the following code to get predictions for the five relevant methods. This runs code, computes node PR, and plots it.

```
python3 main.py -p Wnt -m all --run --node_pr -k 1000 -t 0.5
```

Use `--force` to force the methods to be re-run if necessary.  Final files are here:

```
../plots/node-motivation-BTB-PCSF-PL-RN-RWR-SP-2018-Wnt.pdf
../plots/node-motivation-BTB-PCSF-PL-RN-RWR-SP-2018-Wnt.png
```

## 3. Benchmark figures:

Pick PL as input methods for all pathways; plot composite.

### 3a. DFS/BFS and weighted/unweighted and PerfectLinker Nodes/Edges

```
python3 main.py --benchmark --upper_bounds -p Wnt -m run_PathLinker
python3 main.py --benchmark --upper_bounds -p Wnt -m run_RWR
```

(TODO - there's an error when you run them as one call - must be something with the `--benchmark` option.) Once all the runs and pr curves are calculated, you can re-run the plot with `FULL_AXIS=False` to get a zoomed in view:

```
python3 ../Validation/PR/plot_pr.py /Volumes/compbio/2020-05-PRAUG/runs ../plots PRAUG-PL_2018_Wnt_k500 PRAUG-PL-BFS_2018_Wnt_k500 PRAUG-PL-WEIGHTED_2018_Wnt_k500 PRAUG-PL-BFS-WEIGHTED_2018_Wnt_k500 PRAUG-GT-EDGES_2018_Wnt PRAUG-GT-NODES_2018_Wnt
python3 ../Validation/PR/plot_pr.py /Volumes/compbio/2020-05-PRAUG/runs ../plots PRAUG-RWR_2018_Wnt_k500 PRAUG-RWR-BFS_2018_Wnt_k500 PRAUG-RWR-WEIGHTED_2018_Wnt_k500 PRAUG-RWR-BFS-WEIGHTED_2018_Wnt_k500 PRAUG-GT-EDGES_2018_Wnt PRAUG-GT-NODES_2018_Wnt
```

### 3b. varying k in PL & PRAUG(PL)

Lets set `k=[50,100,500,1000,5000]`
```
python3 main.py -p Wnt --pl_sweep
```

To just plot, set `PARAMS=True` for gray color palette and run:

```
python3 ../Validation/PR/plot_pr.py /Volumes/compbio/2020-05-PRAUG/runs ../plots PRAUG-PL_2018_Wnt_k50 PRAUG-PL_2018_Wnt_k100 PRAUG-PL_2018_Wnt_k500 PRAUG-PL_2018_Wnt_k1000 PRAUG-PL_2018_Wnt_k5000 PL_2018_Wnt_k5000
```

### 3d. varying tau in RWR and PRAUG-RWR.

```
python3 main.py -p Wnt --rwr_sweep
```

To just plot. set `PARAMS=True` for gray color palette and run:

```
python3 ../Validation/PR/plot_pr.py /Volumes/compbio/2020-05-PRAUG/runs ../plots PRAUG-RWR_2018_Wnt_a0.85-t0.1 PRAUG-RWR_2018_Wnt_a0.85-t0.2 PRAUG-RWR_2018_Wnt_a0.85-t0.3 PRAUG-RWR_2018_Wnt_a0.85-t0.4 PRAUG-RWR_2018_Wnt_a0.85-t0.5 PRAUG-RWR_2018_Wnt_a0.85-t0.75 RWR_2018_Wnt_a0.85-t0.75
```

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

## Pathway Case Studies: Notch and Wnt

```
python3 main.py --run --runpraug -p Notch -m all --pr --plot --upper_bounds
python3 main.py --run --runpraug -p Wnt -m all --pr --plot --upper_bounds
```

```
python3 main.py --vote -p Notch --post_graphs aritz@reed.edu platypus
python3 main.py --vote -p Wnt --post_graphs aritz@reed.edu platypus
```
## Others?

### Prediction overlaps (venn or heatmap)
### subset composte on small/med pathways

# Example Runs (for debuging/refactoring)

Anna's runs:
- Reset `DEST_PATH` to be external harddrive.
```
nohup python3 main.py --run --runpraug --pr --plot -p all -m all > /Volumes/compbio/2020-05-PRAUG/runs/2020-05-24.out
 ```
