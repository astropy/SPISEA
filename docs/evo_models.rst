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
* MergedSiessGenevaPadova
* Baraffe15
* Ekstrom12
* Parsec
* Pisa


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

.. autoclass:: evolution.MergedSiessGenevaPadova
	       :show-inheritance:
		  
.. autoclass:: evolution.Baraffe15
	       :show-inheritance:
		  
.. autoclass:: evolution.Ekstrom12
	       :show-inheritance:
		  
.. autoclass:: evolution.Parsec
	       :show-inheritance:	       

.. autoclass:: evolution.Pisa
	       :show-inheritance:
   
