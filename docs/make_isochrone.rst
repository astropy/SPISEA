.. _isochrone_objects:

================
Isochrone Object
================
The isochrone classes are defined in popstar/synthetic.py. The primary
inputs to a isochrone object is the stellar population age, distance,
and metallicity, along with the :ref:`atmo_models` and
:ref:`evo_models`.

If the IsochronePhot sub-class is used, then synthetic photometry
will be produced. Additional inputs such as total extinction,
:ref:`ext_law`, and :ref:`filters` can be defined. 

To define an isochrone::
  
  from popstar import synthetic
  iso = synthetic.<Isochrone_class>

See `Quick Start Example
<https://github.com/astropy/PyPopStar/blob/new_doc/docs/Quick_Start_Make_Cluster.ipynb>`_
for a detailed example for how to make an isochrone and interact with
the output.

Base Isochrone Class
----------------------------
.. autoclass:: synthetic.Isochrone
	       :members: plot_HR_diagram, plot_mass_luminosity



Isochrone Sub-classes
-----------------------

.. autoclass:: synthetic.IsochronePhot
	       :show-inheritance:
		:members: make_photometry, plot_CMD, plot_mass_magnitude
