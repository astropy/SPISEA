.. _evo_models:

========================================
Evolution Model Object
========================================

Stellar evolution models are defined as classes in
popstar/evolution.py. These can be called by::

  from popstar import evolution
  evo = evolution.<model_name>()

The evolution object is an input for the :ref:`isochrone_objects`. 
  
The stellar evolution model classes currently supported by PyPopStar
are:

* MISTv1 (default)
* MergedBaraffePisaEkstromParsec
* MergedPisaEkstromParsec
* Baraffe15
* Ekstrom12
* Parsec
* Pisa

Currently, only the MISTv1 grid offers non-solar metallicity models
(-4 < [Fe/H] <0.5). All other grids have solar metallicity only.
If you require other evolution models or need to
expand the existing grids, please see
:ref:`add_evo_models` and/or let us know on our Github `issue tracker
<https://github.com/astropy/PyPopStar/issues>`_. 


Base Evolution Model Class
------------------------------------
.. autoclass:: evolution.StellarEvolution
	       :members:
 
Specific Evolution Model Classes
--------------------------------------
.. autoclass:: evolution.MISTv1
	       :show-inheritance:


.. autoclass:: evolution.MergedBaraffePisaEkstromParsec
	       :show-inheritance:
		  
.. autoclass:: evolution.MergedPisaEkstromParsec
	       :show-inheritance:	       
		  
.. autoclass:: evolution.Baraffe15
	       :show-inheritance:
		  
.. autoclass:: evolution.Ekstrom12
	       :show-inheritance:
		  
.. autoclass:: evolution.Parsec
	       :show-inheritance:	       

.. autoclass:: evolution.Pisa
	       :show-inheritance:
   
