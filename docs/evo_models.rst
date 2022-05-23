.. _evo_models:

========================================
Evolution Model Object
========================================

Stellar evolution models are defined as classes in
``spisea/evolution.py``. These can be called by::

  from spisea import evolution
  evo = evolution.<model_name>()

The evolution object is an input for the :ref:`isochrone_objects`. 
  
Below is a table of the evolution model grids currently supported by SPISEA.

.. figure:: images/evo_models.png
	    :width: 900
            :height: 172
	    :align: center

Please note the stellar mass, age, and metallicity range of the evolution
model grid you choose. If you require other evolution models or need to
expand the existing grids, please see
:ref:`add_evo_models`. 

ModelMismatch Error
-----------------------------------------
Since SPISEA v2.1.4, we began tracking the version of the evolution
model grid (e.g., the evolution models stored in
``<SPISEA_MODELS>/evolution``).
Each evolution model class has a required model grid
version assigned to it. If your evolution model grid
does not match or exceed the minimum version required by
your desired evolution model, a ``ModelMismatch``
exception will be raised.

To resolve the ``ModelMismatch`` error, please
re-download the latest version of the evolution
model grid in the installation instructions (:ref:`models`).
     

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
		  		  
.. autoclass:: evolution.Baraffe15
	       :show-inheritance:
		  
.. autoclass:: evolution.Ekstrom12
	       :show-inheritance:
		  
.. autoclass:: evolution.Parsec
	       :show-inheritance:	       

.. autoclass:: evolution.Pisa
	       :show-inheritance:
   
