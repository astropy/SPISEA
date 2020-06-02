.. _add_atmo_models:

========================================
Adding Atmosphere Models
========================================
Stellar atmosphere grids are implemented in PyPopStar using STSci's `pysynphot <https://pysynphot.readthedocs.io/en/latest/index.html>`_ and CDBS infrastructure.
Adding a new atmosphere grid involves saving the models in a new sub-directory under your ``PYSYN_CDBS`` directory and creating a new function in ``atmospheres.py`` to access those files.

Detailed documentation on this is coming soon. In the meantime, let us know on the  Github `issue tracker
<https://github.com/astropy/PyPopStar/issues>`_ if you'd like to
implement a new atmospheric model grid.
