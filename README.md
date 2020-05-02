# Pathway Reconstruction

This repository contains the files needed to run and test various methods for
signaling pathway reconstruction, with a focus on HybridLinker.

# Installation Requirements

Most of the code in this repository is written in Python3. A full dependency list
is not available at this time, but at least:

* NetworkX
* Pandas
* Numpy
* GraphSpace

Some of the code in this repository is written in Bash. Because of this, we cannot
guarantee complete functionality on Windows systems at this time.

## Install Python modules

The requirements are listed in `requirements.txt`. To install the python modules to your user directory (in python3), use the command

```
pip3 install --user -r requirements.txt
```

Depending on your install, you might use `pip` instead of `pip3`.

## Install CBC

[CBC](https://github.com/coin-or/Cbc) (Coin-or brank and cut) is an open-source linear programming solver written in C++. It is accessed from Python scripts using the `mip` python module.  CBC is only needed if you plan to run ResponseNet, which relies on a linear program solver. 

On a Mac OS X, the easiest way to install CBC is through Homebrew:

```
brew tap coin-or-tools/coinor
brew install cbc
```

For other platforms, see [their website](https://github.com/coin-or/Cbc) with instructions to download or build from source. 

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
