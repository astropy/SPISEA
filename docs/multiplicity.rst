.. _multi_obj:

===========================
Stellar Multiplicity Object
===========================
The properties of multiple systems in the stellar population is
defined by the stellar multiplicity object. The multiplicity classes
are defined in spisea/imf/multiplicity.py.

To call a multiplicity class::

  from spisea.imf import multiplicity
  multi_obj = multiplicity.<class_name>

The multiplicity object is an input for the :ref:`imf_objects`, as it
impacts how the stellar masses are drawn.  



Stellar Multiplicity Classes
------------------------------------------
.. autoclass:: imf.multiplicity.MultiplicityUnresolved
	       :members: companion_star_fraction,
			 multiplicity_fraction, random_q
