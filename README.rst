====================
PopStar
====================
PopStar is an python package for generating simple stellar populations from synthetic evolution and atmosphere models. Currently, PopStar generates single-age, single-metallicity populations (i.e. star clusters). It supports different initial mass functions, multiplicity perscriptions, reddening laws, filter functions, atmosphere models, and evolution models. The pacakge is object oriented and is extensible. 

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

Path Setup
----------
Add the following to your python path in your .cshrc:

    setenv PYTHONPATH /<path_to_python_direcotries>/PopStar

or, if using bash, add this to your .bashrc:

    export PYTHONPATH=/<path_to_python_direcotries>/PopStar

This path should point to the top-level "PopStar" directory. The actual python functions
should be found in the "popstar" sub-directory immediately below this. 

Install Evolution and Atmosphere Models
---------------------------------------
Once you have cloned the popstar repository, you will need to download the
grid of evolution and atmosphere models. Both sets of models can be found here:

http://w.astro.berkeley.edu/~jlu/popstar/

The required file is poster_models.tar.gz. Note that it is quite large (7.7 GB).

The following evolution models are stored in the "evolution" directory and are ready for 
use:

* baraffe\_pisa\_ekstrom\_parsec (solar metallicity, with and without rotation)
* MIST v1.2 (solar metallicity)

See PopStar/docs/README_EvolutionModels.txt for a description of these models and the
associated references. [under construction]

Evolution models supported by code but not in popstar_models.tar.gz (yet):

* Ekstrom2012 (rot and norot)
* ParsecV1.2S
* Pisa2011
* geneva
* siess\_meynetMaeder_padova
* padova
* pallaStahler1999
* siess2000

The atmosphere models use the CDBS framework (which is supported by
pysynphot). We have added new grids of atmosphere models (including
merged sets of atmospheres). These grids are found in the "cdbs/grid" directory of 
popstar_models.tar.gz. 

Ready: 
* merged\_atlas\_phoenix 
* phoenix\_v16 - high resolution (reformatted for CDBS)
* phoenix\_v16_rebin - downgraded to improve synthetic photometry
performance.

See PopStar/docs/README_AtmoModels.txt for a description of these models sets and the
associated references. [under construction]


Setup Path to Models
--------------------

You need to notify python where these models are going to live. This
is done in two steps.
In your .cshrc or .bashrc, set the PYSYN_CDBS and POPSTAR_MODELS variables:

    setenv PYSYN_CDBS /<path_to_popstar_models_directory>/cdbs
    
    setenv POPSTAR_MODELS /<path_to_popstar_models_directory>/

or

    export PYSYN_CDBS=/<path_to_popstar_models_directory>/cdbs
    
    export POPSTAR_MODELS=/<path_to_popstar_models_directory>

Where the popstar_models directory is wherever you unzipped the popstar_models.tar.gz file. 

Testing PopStar Setup
---------------------
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
    
Otherwise, you should be all set! Try the examples in the Quick Start Guide below to 
make sure everything is working.
    
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
