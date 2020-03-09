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

An example of defining an isochrone::
  
  from popstar import synthetic, evolution
  from popstar import atmospheres, reddening
  import numpy as np

  # Define isochrone parameters
  logAge = np.log10(5*10**6.) # Age in log(years)
  AKs = 0.8 # extinction in Ks-band mags
  dist = 4000 # distance in parsec
  metallicity = 0 # Metallicity in [M/H]

  # Define evolution/atmosphere models and extinction law
  evo_model = evolution.MISTv1() 
  atm_func = atmospheres.get_merged_atmosphere
  red_law = reddening.RedLawHosek18b()

  # Specify filters for synthetic photometry. Here we 
  # use  the HST WFC3-IR F127M, F139M, and F153M filters
  filt_list = ['wfc3,ir,f127m', 'wfc3,ir,f139m', 
  'wfc3,ir,f153m']

  # Make Isochrone object. We will use the IsochronePhot 
  # object since we want synthetic photometry. 
  #
  # Note that is calculation will take a few minutes to run, 
  # unless this isochrone has been generated previously.
  my_iso = synthetic.IsochronePhot(logAge, AKs, dist, 
                            metallicity=0,
                            evo_model=evo_model, 
                            atm_func=atm_func,
                            red_law=red_law, 
                            filters=filt_list)

See `Quick Start Example
<https://github.com/astropy/PyPopStar/blob/new_doc/docs/Quick_Start_Make_Cluster.ipynb>`_
for a detailed example for how to interact with the isochrone object.

Base Isochrone Class
----------------------------
.. autoclass:: synthetic.Isochrone
	       :members: plot_HR_diagram, plot_mass_luminosity



Isochrone Sub-classes
-----------------------

.. autoclass:: synthetic.IsochronePhot
	       :show-inheritance:
		:members: make_photometry, plot_CMD, plot_mass_magnitude
