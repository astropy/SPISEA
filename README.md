# PyPopStar

PyPopStar is an python package for generating simple stellar populations from synthetic evolution and atmosphere models. Currently, PyPopStar generates single-age, single-metallicity populations (i.e. star clusters). It supports different initial mass functions, multiplicity perscriptions, reddening laws, filter functions, atmosphere models, and evolution models. It also produces compact object populations using an initial-final mass relation. The pacakge is object oriented and is extensible. 

One of PyPopStar's strengths is the large set of models that are accessible. In particular, we have created a set of models that "merge" the best-in-class evolution and atmosphere models to cover the full range of masses from 0.08 - 150 Msun at very young ages, including pre-main-sequence evolution.

Here is a brief list of things that PyPopStar can simulation built-in:

* make an isochrone for many different filter sets
* make a star cluster at any age with an unusual IMF
* make a spectrum of a star cluster in integrated light

We encourage contributions to PyPopStar, particular those that add support for star formation histories, new models, higher spectral resolution, etc.


## INSTALL (from git)

### Dependencies

* python (preferably via AstroConda, as it includes some necessary
  packages, like astropy and pysynphot)
* astropy
* pysynphot
* scipy
* numpy
* matplotlib
* STScI CDBS data package (see [here](http://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/synphot-throughput-tables.html) or download directly from ftp://archive.stsci.edu/pub/hst/pysynphot. You need synphot1.tar.gz (63 MB unzipped) and synphot5.tar.gz (791 MB unzipped). See "Set up CDBS directory" below for further instructions.)

### Build

Clone the git repo.
PyPopStar only consists of python codes, so you can either run the
standard `python setup.py install` or you can simply modify your path
environment variables to point to the git python location.

If the installation directory for PyPopStar (where you cloned the repo),
is `<path_to_pypopstar>`, then for BASH, edit your `.bash_profile` to
include:

    export PYTHONPATH=$PYTHONPATH:<path_to_pypopstar>/PyPopStar

where you replace `<path_to_pypostar>` with the appropriate directory. 

For C-shell, edit your .cshrc to include:

    setenv PYTHONPATH <path_to_pypopstar>/PyPopStar

where you replace `<path_to_postar>` with the appropriate
directory. Note, the python path should point to the top-level "PyPopStar"
directory. The actual python functions should be found in the
`PyPopStar/popstar/` sub-directory.

### Set up CDBS directory
PyPopStar uses the STScI CDBS infrastructure to store model atmospheres and HST filter functions. It expects a top-level `cdbs` directory with 3 sub-directories: `comp`, `grid`, and `mtab`. The `comp` and `mtab` directories are obtained from the `synphot1.tar.gz` file you downloaded (grp/hst/cdbs/comp, grp/hst/cdbs/mtab). The `grid` directory will eventually hold all of the model atmospheres, but start by putting the Phoenix model directory from `synphot5.tar.gz` there (grp/hst/cdbs/grid/phoenix). You will add additional atmosphere model directories to `cdbs/grid` in the "Install Evolution and Atmosphere Models" section below.

So, you should end up with a single `cdbs` directory with the `cdbs/comp`, `cdbs/grid`, and `cdbs/mtab` directories defined.

### Install Evolution and Atmosphere Models

Once you have cloned the PyPopStar repository, you will need to download the
grids of evolution and atmosphere models. The evolution models are
stand-alone and specific to PyPopStar. The atmosphere models use the
STScI CDBS conventions and should be copied into your `cdbs/grid` directory.

The two files you need to download are:

[popstar_models.tar.gz](http://astro.berkeley.edu/~jlu/popstar/popstar_models.tar.gz)  (2.6 GB; 18 GB unzipped)

[postar_cdbs.tar.gz](http://astro.berkeley.edu/~jlu/popstar/popstar_cdbs.tar.gz)  (142 MB; 248 MB unzipped)

You can also optionally download a third file, which contains high-resolution versions of the atmospheres in `popstar_cdbs.tar.gz`:

[popstar_cdbs_highres.tar.gz](http://astro.berkeley.edu/~jlu/popstar/popstar_cdbs_highres.tar.gz) (50 GB; 74 GB unzipped)

However, `popstar_cdbs_highres.tar.gz` is not needed for PyPopStar to run. PyPopStar uses the low-resolution atmospheres in `popstar_cdbs.tar.gz` for synthetic photometry by default, as this is much faster and is sufficient in most applications. 

Once downloaded, `popstar_cdbs.tar.gz` (and `popstar_cdbs_highres.tar.gz`, if desired) should be expanded in 
the directory that houses `cdbs`. The `popstar_models.tar.gz` file can be expanded anywhere, as you will set
an environment variable to its location. However, we recommend downloading it to the same location as your cdbs directory. 

For instance, if your cdbs directory is `/g/lu/models/cdbs/` then you should put the .gz files in `/g/lu/models/` and then unzip them from there: 


```console
cd /g/lu/models/
tar xvf popstar_cdbs.tar.gz
tar xvf popstar_models.tar.gz
```

`popstar_cdbs.tar.gz` will put the model atmospheres in `cdbs/grid`, and `popstar_models.tar.gz` will put the evolution models in a directory called `evolution/`. 

See [README_EvolutionModels.txt](https://github.com/astropy/PyPopStar/blob/master/docs/README_EvolutionModels.txt) for a description of the stellar evolution models and [README_AtmoModels.txt](https://github.com/astropy/PyPopStar/blob/master/docs/README_AtmoModels.txt) for a description of the stellar atmosphere models included in this download.
You can add additional evolution models and atmosphere models to PyPopStar if you wish.

### Setup Path to Models

You need to notify python where the evolution and atmosphere models live. This
is done with two environment variables, PYSYN_CDBS and POPSTAR_MODELS. 
In your .cshrc or .bashrc, set them as the following:

```sh
setenv PYSYN_CDBS /<path_to_models_directory>/cdbs
setenv POPSTAR_MODELS /<path_models_directory>/
```

or

```sh
export PYSYN_CDBS=/<path_to_models_directory>/cdbs
export POPSTAR_MODELS=/<path_to_models_directory>
```

Where the <path_to_models_directory> is where you downloaded and expanded the `popstar_cdbs.tar.gz` and `popstar_models.tar.gz` files (e.g., where the `cdbs/` and `evolution/` directories live). 

## Testing PyPopStar Setup

If all goes well, you should be able to import any of the PyPopStar functions
in a python environment window using an import statement like those at the top
of the Quick Start Guide, such as:
    
    from popstar import synthetic
    
If the `PYTHONPATH` is broken, then python won't know where to find the popstar codes and
you will get an error reporting "No module named popstar". If the `POPSTAR_MODELS` or 
`PYSYN_CDBS` paths are broken, then pypopstar won't know here to go to get the 
stellar models. and you will get either or both of the following warnings:

    UserWarning: PYSYN_CDBS is undefined; functionality will be SEVERELY crippled.
    
    UserWarning: POPSTAR_MODELS is undefined; functionality will be SEVERELY crippled.
    
Otherwise, you should be all set! You may get warnings that say:

    UserWarning: Extinction files not found in /u/mwhosek/models/cdbs/extinction
    
    UserWarning: No thermal tables found, no thermal calculations can be performed
    
But these can be safely ignored since PyPopstar doesn't use this functionality.

To further test the installation, try running the Quick Start Guide [here](https://github.com/astropy/PyPopStar/blob/master/docs/Quick_Start_Make_Cluster.ipynb). 
To test the range of evolution models, atmosphere models, and photometric
filters available in PyPopStar, run the `test_evolution_models()`, `test_atmospheres_models()`, and `test_filters()` functions in `popstar/tests/test_models.py`. 
    
## Documentation

For a quick tutorial on how to make a star cluster with pypopstar, see
the jupyter notebook at `PyPopstar/docs/Quick_Start_Make_Cluster.ipynb`.
Additional Jupyter notebooks with tutorials to produce the plots shown in the PyPopStar paper (Hosek et al., in prep)
can be found in `PyPopStar/docs/paper_examples/`. 

Additional documentation:

* Stellar Evolution Models: [docs/README_EvolutionModels.txt](https://github.com/astropy/PyPopStar/blob/master/docs/README_EvolutionModels.txt)
* Stellar Atmosphere Models: [docs/README_AtmoModels.txt](https://github.com/astropy/PyPopStar/blob/master/docs/README_AtmoModels.txt)
* Extinction: [docs/README_Extinction.txt](https://github.com/astropy/PyPopStar/blob/master/docs/README_Extinction.txt)
* Filters: [docs/README_Filters.txt](https://github.com/astropy/PyPopStar/blob/master/docs/README_Filters.txt)
* Initial Mass Function: docs/README_IMF.txt [under construction]
* Multiplicity: docs/README_Multiplicity.txt [under construction]
* Initial-Final Mass Relation: docs/README_IFMR.txt [under construction]


## Other Resources

* _Astropy: http://www.astropy.org/
* _git: http://git-scm.com/
* _github: http://github.com
* _Cython: http://cython.org/

## License 
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
