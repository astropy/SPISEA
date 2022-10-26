.. sampledoc documentation main file, created by
   sphinx-quickstart on Tue Aug 11 05:04:40 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. figure:: ../SPISEA_logo_final.png
	    :width: 380
            :height: 336
	    :align: center

==========================================================================
SPISEA: Stellar Population Interface for Stellar Evolution and Atmospheres
==========================================================================
SPISEA is an open-source python package that generates single-age, single-metallicity
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

SPISEA can be downloaded from `Github <https://github.com/astropy/SPISEA>`_. Please cite `Hosek et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200606691H/abstract>`_  if you use SPISEA in
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
Have an idea? We'd love to have you join our team!
We welcome and encourage contributions to SPISEA.

For feature additions, we ask that users make their
own fork of the repository, make their changes, and then submit a pull
request to the "dev" branch.

All contributions will be acknowledged on the :ref:`contributors` page
(with permission). Contributors with features used in code
releases will be co-authors in future SPISEA software papers.


Change Log
----------

2.1.7 (2022-10-26)
  * Made several updates to the extinction law objects, including:
    
    * Added Indebetouw05, Schoedel10, Fritz11, and NoguerasLara20 laws
    * Improved docs and added plotting functions for Nishiyama09 and RomanZuniga07 laws
    * Updated K-band scale wavelength for Fitzpatrick07 and NoguerasLara18 to match published values
    * Updated assorted "under the hood" functions to handle new laws
    * Added more test functions for the extinction laws
    * Updated :ref:`ext_law` page so laws are in chronological order

2.1.6 (2022-09-15)
  * Bugfix in how evolution object identifies what age and metallicity
    model to pull from the grid of available models. Previously, it used
    first grid model where grid model value >
    requested value (numpy.searchsorted routine). Now, is takes grid model
    where the difference between the grid model value and requested
    value is minimized. This is most relevant for metallicity
    selection for MIST isochrones, where the model grid metallicities are
    spaced out. 

2.1.5 (2022-07-20)
  * Added additional parameters to get_bb_atmosphere function to give
    user control over the min and max wavelength values, as well as
    wavelength sampling. See :ref:`atmo_models` for details.

2.1.4 (2022-05-24)

  * Expanded MISTv1 evolution model grid age range to logAge = 5.01 --
    10.30
  * Begun version tracking for evolution model grid; updated grid
    (with new MISTv1 age range) is grid 1.1
  * A new ModelMismatch
    exception will be raised if evolution model object requires a
    higher version of the model grid than the user has installed. See
    :ref:`evo_models` for details
  * Fixed MISTv1 sub-directory names so proper Z_solar value
    of 0.0142 is used (only naming convention, does not effect output
    results)
  * Updated evolution model and atmosphere model doc pages with tables
    describing properties of installed options
  * Small tweaks to installation instructions for clarity
  * Updated test functions and removed outdated dependence on nose package
 
2.1.3 (2022-04-23)

  * Updates to documentation and Quick-Start Guide
  * Fully depreciated Hosek18 extinction law to avoid confusion with Hosek18b

2.1.2 (2022-04-01)

  * New installation methods using Docker containers and Singularity


2.1.1 (2022-03-30)

  * Bugfix: fixed application of IFMR to objects in companions table
    (ResolvedCluster and associated subclasses: applies when
    multiplicity + IFMR are used)

  * Minor additions to docs

2.1.0 (2022-01-24)

  * Added ability to create resolved multiple systems, with orbital
    elements assigned to each companion (see :ref:`multi_obj`)

  * Added new metallicity-dependent IFMR models (see
    :ref:`ifmr_objects`)

  * Added blackbody atmosphere class

  * Additional test functions, filter functions, and documentation
    updates

2.0.0 (2020-07-09)

  * Top-level software name change from PyPopStar to SPISEA (see
    :ref:`version` for instructions to do this). This version of code was released in `Hosek et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200606691H/abstract>`_. 

1.0.1 (2020-06-24)

  * Bug fix for photometric column headers for some filters, added new tests regarding total cluster mass, small documentation edits

1.0.0 (2019-12-01)

  * Initial release


Logo by Natasha Abrams
