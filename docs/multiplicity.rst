.. _multi_obj:

===========================
Stellar Multiplicity Object
===========================
The properties of multiple systems in the stellar population is
defined by the stellar multiplicity object. The multiplicity classes
are defined in ``spisea/imf/multiplicity.py``.

To call a multiplicity class::

  from spisea.imf import multiplicity
  multi_obj = multiplicity.<class_name>

The multiplicity object is an input for the :ref:`imf_objects`, as it
impacts how the stellar masses are drawn.

The user can choose either
an unresolved or a resolved multiplicity object. If a resolved
multiplicity object is selected, then orbital parameters are
assigned to each companion star (e.g semi-major axis, eccentricity,
inclination). These values are added as additional columns in the ``companions``
table off of the cluster object. Note that the synthetic photometry
for the system is unchanged, i.e., the photometry reported in the
``star_systems`` table off the cluster object still reports the combined photometry of the entire
system. 



Unresolved Multiplicity Classes
------------------------------------------
.. autoclass:: imf.multiplicity.MultiplicityUnresolved
	       :members: companion_star_fraction,
			 multiplicity_fraction, random_q


Resolved Multiplicity Classes
------------------------------------------
.. autoclass:: imf.multiplicity.MultiplicityResolvedDK
	       :show-inheritance:
