.. _add_evo_models:

========================================
How to Add Evolution Models
========================================
Need an evolution model not listed among those pre-packaged
with PyPopStar? No problem! With the additional documentation provided
below, you can implement your own evolution models into the PyPopStar
framework. 

The PyPopStar file containing the code associated with the evolution
models is ``evolution.py``. Each evolution model grid is implemented
as a sub-class of the base ``StellarEvolution`` object: the ``MISTv1``
sub-class corresponds to the MISTv1 model grid, the ``Parsec`` sub-class is
for the Parsec model grid, etc. If you are planning to add to an
additional model grid, you will need to edit the appropriate
sub-class. If you want to add a brand new model grid, you will need to
create a new sub-class for that model.

If you have questions or run into problems, please raise an issue on
our Github `issue tracker
<https://github.com/astropy/PyPopStar/issues>`_. If you are willing to
have the new models you add be added to the PyPopStar package and made
available to the community, please fork or branch off of the
development repository and then submit merge requests to add your
changes to the package. These contributions will be added to the
master branch (with attributes given!) in the next code update.


Adding New Metallicities to An Existing Model Grid
--------------------------------------------------








Adding New Ages to An Existing Model Grid
--------------------------------------------------
