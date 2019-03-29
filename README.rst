====================
PopStar
====================
PopStar is an astropy affiliated package for generating simple stellar populations from synthetic evolution and atmosphere models. Currently, PopStar generates single-age, single-metallicity populations (i.e. star clusters). It supports different initial mass functions, multiplicity perscriptions, reddening laws, filter functions, atmosphere models, and evolution models. The pacakge is object oriented and is extensible. 

PopStar's primary strength is the large set of models that are accessible. In particular, we have created a set of models that "merge" the best-in-class evolution and atmosphere models to cover the full range of masses from 0.08 - 150 Msun at very young ages, including pre-main-sequence evolution.

Here is a brief list of things that PopStar can simulation built-in:

* make an isochrone for many different filter sets
* make a star cluster at any age with an unusual IMF
* make a spectrum of a star cluster in integrated light

We encourage contributions to PopStar, particular those that add support for star formation histories, metallicity distributions, new models, higher spectral resolution, etc.


INSTALL (from git)
==================

Dependencies
------------
* python (preferably via AstroConda, as it includes some necessary
  packages, like astropy and pysynphot)
* astropy
* pysynphot
* scipy
* numpy
* matplotlib
* STScI CDBS data package

Build
----------
PopStar only consists of python codes, so you can either run the
standard `python setup.py install` or you can simply modify your path
environment variables to point to the git python location.

If the installation directory for PopStar (where you cloned the repo),
is `<path_to_popstar>`, then for BASH, edit your `.bash_profile` to include:
```sh
export PYTHONPATH=$PYTHONPATH:<path_to_popstar>/PopStar
```
where you replace `<path_to_postar>` with the appropriate directory. 

For C-shell, edit your `.cshrc` to include:

```sh
setenv PYTHONPATH <path_to_popstar>/PopStar
```
where you replace `<path_to_postar>` with the appropriate directory. 


Install Evolution and Atmosphere Models
---------------------------------------
Once you have cloned the popstar repository, you will need to download the
grid of evolution and atmosphere models. The evolution models are
stand-alone and specific to PopStar. The atmosphere models use the
STScI CDBS conventions and should be stored in your already installed
CDBS directory.

The two files to download (but not yet expand) are:
<http://astro.berkeley.edu/~jlu/popstar/popstar_models.tar.gz>
<http://astro.berkeley.edu/~jlu/popstar/popstar_models.tar.gz>

The `popstar_cdbs.tar` file should be expanded in the directory that
houses `cdbs`.
The `popstar_models.tar` file can be expanded anywhere as you will set
an environment variable to the install location. However, we recommend
installing it in the same location as your cdbs. 
For instance, if your cdbs install is in
`/g/lu/models/cdbs/` then you should:

```sh
cd /g/lu/models/
tar xvf popstar_cdbs.tar
tar xvf popstar_models.tar
```

Next, you will need to set environment variabls to point to these
model directories. See Setup Path to Models below.


Setup Path to Models
--------------------

You need to notify python where these models are going to live. This
is done in two steps.

In your .cshrc or .bashrc, set the PYSYN_CDBS and POPSTAR_MODELS variables:
```sh
setenv PYSYN_CDBS /g/lu/models/cdbs
setenv POPSTAR_MODELS /g/lu/models
```
or
```sh
export PYSYN_CDBS=/g/lu/models/cdbs
export POPSTAR_MODELS=/g/lu/models
```

Quick Start Guide
-------------------
For a quick tutorial on how to make a star cluster with popstar, see
the jupyter notebook at Popstar/docs/Quick_Start_Make_Cluster.ipynb
    

Other Resources
===============

* _Astropy: http://www.astropy.org/
* _git: http://git-scm.com/
* _github: http://github.com
* _Cython: http://cython.org/
