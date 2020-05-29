.. _add_evo_models:

========================================
How to Add Evolution Models
========================================
Need an evolution model not listed among those pre-packaged
with PyPopStar? No problem! With the additional documentation provided
below, you can implement your own evolution models into the PyPopStar
framework.

If you have questions or run into problems, please raise an issue on
our Github `issue tracker
<https://github.com/astropy/PyPopStar/issues>`_. If you are willing to
have the new models you add be added to the PyPopStar package and made
available to the community, please fork or branch off of the
development repository and then submit merge requests to add your
changes to the package. These contributions will be added to the
master branch (with attributes given!) in the next code update.

Setup: Codes, Directories, Files
--------------------------------
The codes associated with the evolution
models are in ``evolution.py``. Each evolution model grid is implemented
as a sub-class of the ``StellarEvolution`` object: the ``MISTv1``
sub-class corresponds to the MISTv1 model grid, the ``Parsec`` sub-class is
for the Parsec model grid, etc. If you are planning to add to an
additional model grid, you will need to edit the existing
sub-class. If you want to add a brand new model grid, you will need to
create a new sub-class for that model.

The evolution model grids are packaged as individual isochrone files
in PyPopStar. The original source of these files are typically the
website for the given evolution model, which is provided in the
documentation of the different evolution model sub-classes.
The isochrone files are stored in the
``$POPSTAR_MODELS/evolution/`` directory. Each evolution model grid is contained
in its own sub-directory; for example,
``$POPSTAR_MODELS/evolution/PARSECV1.2s`` contains the PARSEC
isochrones.
Each of these sub-directories has a similar structure::

  ---> version (if applicable)
              ---> iso
	                  ---> metallicity
			              ---> rotation (if applicable)
				                  ---> individual files

The metallicity directories are named ``z<mass_fraction>``, where
``<mass_fraction>`` is the mass fraction of metals. The individual
isochrone files saved as FITS tables with the name convention
``iso_<logAge>fits``, where ``<logAge>`` is the log(age) of the
population (age in years). While the ``<mass_fraction>`` can be any
number of digits you want, PyPopStar expects the ``<logAge>`` to
always have two places after the decimal.




Adding New Metallicities to An Existing Model Grid
--------------------------------------------------








Adding New Ages to An Existing Model Grid
--------------------------------------------------



Creating a New Model Grid
-------------------------
