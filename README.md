# Pathway Reconstruction

This repository contains the files needed to run and test various methods for
signaling pathway reconstruction, with a focus on HybridLinker.

Most of the code in this repository is written in Python3. A full dependency list
is not available at this time, but at least:

* NetworkX
* Pandas
* Numpy
* GraphSpace

Some of the code in this repository is written in Bash. Because of this, we cannot
guarantee complete functionality on Windows systems at this time.


# Instructions for checking out this repo.

There are three submodules in this repo.  After cloning, you need to update the submodules with

```
git submodule init
git submodule update
```

Check the `.gitmodules` file for more details.

## Adding the interactome.

Within the `localised-pathlinker` directory, the interactome is in a zipped file. We need to unzip it and move it to a place where all methods expect to see it.

Locate and unzip the file.
```
cd Reconstruction/Methods/localized-pathlinker/Data
unzip PathLinker_2018_human-ppi-weighted-cap0_75.txt.zip
```

Make an `Interactomes/` directory and move the interactome. Unfortunately a symbolic link won't work here.
```
mkdir ../../../../Interactomes/
mv PathLinker_2018_human-ppi-weighted-cap0_75.txt ../../../../Interactomes/
```
