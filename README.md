# PopStar

PopStar is an python package for generating simple stellar populations from synthetic evolution and atmosphere models. Currently, PopStar generates single-age, single-metallicity populations (i.e. star clusters). It supports different initial mass functions, multiplicity perscriptions, reddening laws, filter functions, atmosphere models, and evolution models. The pacakge is object oriented and is extensible. 

PopStar's primary strength is the large set of models that are accessible. In particular, we have created a set of models that "merge" the best-in-class evolution and atmosphere models to cover the full range of masses from 0.08 - 150 Msun at very young ages, including pre-main-sequence evolution.

Here is a brief list of things that PopStar can simulation built-in:

* make an isochrone for many different filter sets
* make a star cluster at any age with an unusual IMF
* make a spectrum of a star cluster in integrated light

We encourage contributions to PopStar, particular those that add support for star formation histories, metallicity distributions, new models, higher spectral resolution, etc.


## INSTALL (from git)

### Dependencies

* python (preferably via AstroConda, as it includes some necessary
  packages, like astropy and pysynphot)
* astropy
* pysynphot
* scipy
* numpy
* matplotlib
* STScI CDBS data package (see [here](http://www.stsci.edu/hst/observatory/crds/throughput.html); you need cdbs/comp and cdbs/mtab sub directories (synphot1.tar.gz) as well as phoenix model atmospheres (synphot5.tar.gz) located in the cdbs/grid subdirectory.)

### Build

Clone the git repo.
PopStar only consists of python codes, so you can either run the
standard `python setup.py install` or you can simply modify your path
environment variables to point to the git python location.

If the installation directory for PopStar (where you cloned the repo),
is `<path_to_popstar>`, then for BASH, edit your `.bash_profile` to
include:

    export PYTHONPATH=$PYTHONPATH:<path_to_popstar>/PopStar

where you replace `<path_to_postar>` with the appropriate directory. 

For C-shell, edit your .cshrc to include:

    setenv PYTHONPATH <path_to_popstar>/PopStar

where you replace `<path_to_postar>` with the appropriate
directory. Note, the python path should point to the top-level "PopStar"
directory. The actual python functions should be found in the
`PopStar/popstar/` sub-directory.


### Install Evolution and Atmosphere Models

Once you have cloned the popstar repository, you will need to download the
grid of evolution and atmosphere models. The evolution models are
stand-alone and specific to PopStar. The atmosphere models use the
STScI CDBS conventions and should be copied into your already installed
CDBS directory.

The two files to download (but not yet expand) are:

[popstar_models.tar.gz](http://astro.berkeley.edu/~jlu/popstar/popstar_models.tar.gz)  (2.6 GB)

[postar_cdbs.tar.gz](http://astro.berkeley.edu/~jlu/popstar/popstar_cdbs.tar.gz)  (7 GB)

The `popstar_cdbs.tar` file should be expanded in the directory that
houses `cdbs`.
The `popstar_models.tar` file can be expanded anywhere as you will set
an environment variable to the install location. However, we recommend
installing it in the same location as your cdbs. 
For instance, if your cdbs install is in
`/g/lu/models/cdbs/` then you should:


```console
cd /g/lu/models/
tar xvf popstar_cdbs.tar
tar xvf popstar_models.tar
```

See PopStar/docs/README_EvolutionModels.txt for a description of these
models and the associated references, as well as other evolution
models supported by the code.

The list of supported evolution models includes

* merged grid: ATLAS9, PHOENIXv16, BTSettl
* ATLAS9: Castelli & Kurucz 2004
* PHOENIXv16: Husser et al. 2013
* BTSettl: Allard et al.
* Koester et al. 2010 (white dwarfs)
* Kurucz 1993
  
See PopStar/docs/README_AtmoModels.txt for a description of these
models sets and the associated references. 


### Setup Path to Models


You need to notify python where these models are going to live. This
is done in two steps.
In your .cshrc or .bashrc, set the PYSYN_CDBS and POPSTAR_MODELS variables:

```sh
setenv PYSYN_CDBS /<path_to_models_directory>/cdbs
setenv POPSTAR_MODELS /<path_models_directory>/
```

or

```sh
export PYSYN_CDBS=/<path_to_models_directory>/cdbs
export POPSTAR_MODELS=/<path_to_models_directory>
```

Where the models directory contains the PopStar `evolution/` and CDBS
`cdbs/grid` directories with PopStar atmospheres in it. Your CDBS should
have 3 sub-directories: cdbs/grid, cdbs/comp, and cdbs/mtab. The comp and 
mtab directories come from the CDBS install. 

## Testing PopStar Setup

If all goes well, you should be able to import any of the PopStar functions
in a python environment window using an import statement like those at the top
of the Quick Start Guide, e.g.:
    
    from popstar import synthetic
    
If the PYTHONPATH is broken, then python won't know where to find the popstar codes and
you will get an error reporting "No module named popstar". If the POPSTAR_MODELS or 
PYSYN_CDBS paths are broken, then popstar won't know here to go to get the 
stellar models. In this case, upon import you will get either or both of 
the following warnings:

    UserWarning: PYSYN_CDBS is undefined; functionality will be SEVERELY crippled.
    
    UserWarning: POPSTAR_MODELS is undefined; functionality will be SEVERELY crippled.
    
Otherwise, you should be all set! You may get warnings that say:

    UserWarning: Extinction files not found in /u/mwhosek/models/cdbs/extinction
    
    UserWarning: No thermal tables found, no thermal calculations can be performed
    
But these can be safely ignored. Try the examples in the Quick Start Guide below to 
make sure everything is working.
    
## Documentation

For a quick tutorial on how to make a star cluster with popstar, see
the jupyter notebook at Popstar/docs/Quick_Start_Make_Cluster.ipynb.

Additional documentation:

* Stellar Evolution Models: docs/README_EvolutionModels.txt
* Stellar Atmosphere Models: docs/README_AtmoModels.txt
* Extinction: docs/README_Extinction.txt 
* Filters: docs/README_Filters.txt 
* Initial Mass Function: docs/README_IMF.txt [under construction]
* Multiplicity: docs/README_Multiplicity.txt [under construction]
* Initial-Final Mass Relation: docs/README_IFMR.txt [under construction]


## Other Resources

* _Astropy: http://www.astropy.org/
* _git: http://git-scm.com/
* _github: http://github.com
* _Cython: http://cython.org/
