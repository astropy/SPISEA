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


## INSTALL (from git)

### Dependencies
* python (preferably via AstroConda, as it includes some necessary
  packages, like astropy and pysynphot)
* astropy
* pysynphot
* scipy
* numpy
* matplotlib

### Path Setup      
Add all of these to your python path in your .cshrc:

    set CODE_DIR = /path/to/top/level/python/code
    setenv PYTHONPATH ${CODE_DIR}/popstar/PopStar

### Install Evolution and Atmosphere Models
Once you have cloned the popstar repository, you will need to download the
grid of evolution and atmosphere models. The evolution models can be
found here:

GIT SITE -- load these up
/g/lu/models/evolution/

Ready:

* pisa\_ekstrom_parsec
* pisa\_ekstrom_parsec/norot

Data files ready but not yet supported in code:

* Ekstrom2012 (rot and norot)
* ParsecV1.2S
* Pisa2011
* geneva
* merged
* siess\_meynetMaeder_padova
* padova
* pallaStahler1999
* siess2000

Update:

* STERN2011

The atmosphere models use the CDBS framework (which is supported by
pysynphot). We have added new grids of atmosphere models (including
merged sets of atmospheres). You can find the complete CDBS grid you
can download from here:

GIT SITE -- load these up
/g/lu/models/cdbs

Ready: 

* phoenix\_v16 - high resolution (reformatted for CDBS)
* phoenix\_v16_rebin - downgraded to improve synthetic photometry
performance.


#### Setup Path to Models

You need to notify python where these models are going to live. This
is done in two steps.

In your .cshrc or .bashrc, set the PYSYN_CDBS and POPSTAR_MODELS variables:

    setenv PYSYN_CDBS /g/lu/models/cdbs
    setenv POPSTAR_MODELS /g/lu/models

or

    export PYSYN_CDBS=/g/lu/models/cdbs
    export POPSTAR_MODELS=/g/lu/models


## Other Resources

* _Astropy: http://www.astropy.org/
* _git: http://git-scm.com/
* _github: http://github.com
* _Cython: http://cython.org/
