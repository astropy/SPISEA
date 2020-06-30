.. sampledoc documentation master file, created by
   sphinx-quickstart on Tue Aug 11 05:04:40 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SPISEA (Star Population Interface for Stellar Evolution and
Atmospheres) is an python package that generates single-age, single-metallicity
populations (i.e. star clusters). It gives the user control over many parameters:

* Cluster characteristics (age, metallicity, mass, distance)
* Total extinction, differential extinction, and extinction law
* Stellar evolution and atmosphere models
* Stellar multiplicity and Initial Mass Function
* Initial-Final Mass Relation
* Photometric filters

Here is a brief list of things that SPISEA can do:

* make a cluster isochrone in many filters using different stellar models
  
* make a star cluster at any age with an unusual IMF and unresolved multiplicity
  
* make a spectrum of a star cluster in integrated light

Please cite `Hosek et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200606691H/abstract>`_  if you use SPISEA in
your research.

Getting Started
----------------
.. toctree::
   :maxdepth: 1

   getting_started.rst
   quick_start.rst
   more_examples.rst


Documentation
-------------------------
.. toctree::
   :maxdepth: 1

   make_isochrone.rst
   make_cluster.rst
   evo_models.rst
   atmo_models.rst
   imf.rst
   multiplicity.rst
   ifmr.rst
   extinction.rst
   filters.rst

Advanced Documentation
-----------------------------
.. toctree::
   :maxdepth: 1
	      
   add_evo_models.rst
   add_atmo_models.rst
   add_filters.rst

Contributions
---------------
We encourage contributions to SPISEA, particular
those that add support for star formation histories,
new models, higher spectral resolution, etc.
For feature additions, we ask that users fork or
branch off of the development repository, make their changes,
and then submit merge and pull requests.
