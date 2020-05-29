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

Setup: Codes and Files
----------------------
The codes associated with the evolution
models are in ``evolution.py``. Each evolution model grid is implemented
as a sub-class of the ``StellarEvolution`` object: the ``MISTv1``
sub-class corresponds to the MISTv1 model grid, the ``Parsec`` sub-class is
for the Parsec model grid, etc. If you are planning to add to an
additional model grid, you will need to edit the existing
sub-class. If you want to add a brand new model grid, you will need to
create a new sub-class for that model.

The files containing the evolution model grids are stored in the
``$POPSTAR_MODELS/evolution/`` directory. Each model set is contained
in its own sub-directory; for example, ``$POPSTAR_MODELS/evolution/PARSECV1.2s`` contains the PARSEC evolution
models. Each of these sub-directories has a similar structure:








Adding New Metallicities to An Existing Model Grid
--------------------------------------------------








Adding New Ages to An Existing Model Grid
--------------------------------------------------



Creating a New Model Grid
-------------------------
