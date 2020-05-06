# Pathway Reconstruction

This repository contains the files needed to run and test various methods for
signaling pathway reconstruction, with a focus on HybridLinker. There are thee things you need to do to get this code up and running:

1. Update the repo submodules
2. Put the interactome in the proper place
3. Install dependencies (including method-specific dependencies if you plan to run them)

# Instructions for checking out this repo

There are three submodules in this repo.  After cloning, you need to update the submodules with

```
git submodule init
git submodule update
```

Check the `.gitmodules` file for more details.  

# Adding the interactome

Within the `localized-pathlinker` directory, the interactome is in a zipped file. We need to unzip it and move it to a place where all methods expect to see it.

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

Alternatively, install dependencies in a new environment using `conda` with

```
conda env create -f depends.yml
```

## Method-Specific Installation Requirements


### Install CBC (ResponseNet)

[CBC](https://github.com/coin-or/Cbc) (Coin-or brank and cut) is an open-source linear programming solver written in C++. It is accessed from Python scripts using the `mip` python module.  CBC is only needed if you plan to run ResponseNet, which relies on a linear program solver. 

On a Mac OS X, the easiest way to install CBC is through Homebrew:

```
brew tap coin-or-tools/coinor
brew install cbc
```

For other platforms, see [their website](https://github.com/coin-or/Cbc) with instructions to download or build from source. 

### Install the [Omics Integrator 2](https://github.com/fraenkel-lab/OmicsIntegrator2) (PCSF)

The most recent implementation of the Prize Collecting Steiner Forest (PCSF) is within the `Forest` module of the Omics Integrator v2.  This is the maintained version of the code, which can be installed by cloning the repo and running `setup.py`. This is **not** included as a submodule for the reasons listed below; so make sure the `OmicsIntegrator2` is included in your `PYTHONPATH`.

```
git clone git@github.com:fraenkel-lab/OmicsIntegrator2.git
cd OmicsIntegrator2/
python3 setup.py install --user
```

Depending on your distribution, this may not work the first time.  We discovered that there is no `pcst_fast` module in Python3 (though it exists in Python2 via pypi).  To remedy this, we 

1. **Installed `pcst_fast` from source**, fixing the install so Python3 libraries are called.

```
git clone git@github.com:fraenkel-lab/pcst_fast.git
```

In the `Makefile`, I had to make the following changes to the `pcst_fast_py` line so python3 libraries were called:
- Added `-undefined dynamic_lookup`
- Changed `python-config` to `python3-config`

The line now looks like:
```
pcst_fast_py: $(PCST_FAST_PY_SRC_DEPS:%=$(SRCDIR)/%)
	$(CXX) $(CXXFLAGS) -undefined dynamic_lookup -shared -I $(SRCDIR) -I external/pybind11/include `python3-config --cflags --ldflags` $(SRCDIR)/pcst_fast_pybind.cc $(SRCDIR)/pcst_fast.cc -o pcst_fast.so
```

You can then install `pcst_fast` with

```
make pcst_fast_py
```

You should now be able to `import pcst_fast` in your Python environment. 

2. **Added `pcst_fast` to the `PYTHONPATH` environment variable.** In your `~/.bashrc` or `~/.bash_profile` file (for linux and macs, respectively):

``
export PYTHONPATH=$PYTHONPATH:/Users/aritz/Documents/github/pcst_fast/
```

Remember to `source ~/.bashrc` or `source ~/.bash_profile` for the change to take effect.

3. **Disable `pcst_fast` from the requirements in `OmicsIntegrator`** (since we satisfied the requirement with #1 and #2). In the `setup.py` file, comment out the `pcst_fast` requirement. Lines 19-29 should now look like

```
install_requires=[
        "numpy",
        "pandas==0.23.4",
        "networkx==2.1",
        #"pcst_fast==1.0.7",
        "python-louvain",
        "goenrich",
        "sklearn",
        "axial",
        "scipy"
    ],
```


In the `OmicsIntegrator/` directory, run

```
python3 setup.py install --user
```

(Alternatively, you could comment the `pcst_fast` requirement in `requirements.txt` and install those dependencies. We didn't check this).  

You should now be able to `import OmicsIntegrator` in your Python environment.


