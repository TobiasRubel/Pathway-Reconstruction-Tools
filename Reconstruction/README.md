# Pipelines for Generating Figures

## Generating the composite figures

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

## Generating the node motivation figure

1. Run the following code to get predictions for the five relevant methods. (TODO need to add RWR)

```
python3 main.py -p Wnt -m run_PCSF run_ResponseNet run_BowtieBuilder run_ShortestPaths run_PathLinker
```

2. Run the node PR for these, subsampling 50X nodes and also ignoring adjacents: (TODO need to add RN, PCSF)

```
python3 main.py --node_pr -p Wnt
```

3. Plot the methods.

```
python3 ../Validation/PR/plot_pr.py ../Validation/PR/data ../../../plots ShortestPaths_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000 PathLinker_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000 BowtieBuilder_PathLinker_2018_human-ppi-weighted-cap0_75_Wnt_10000
```

Final file is here:

```
../../../plots/node-motivation-BowtieBuilder-PathLinker-ShortestPaths-PathLinker-2018-human-ppi-weighted-cap0-75-Wnt-10000.png
```
