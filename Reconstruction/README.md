## Generating the node motivation figure.

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
