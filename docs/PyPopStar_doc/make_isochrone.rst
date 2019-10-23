.. _isochrone_objects:

******************************************
Isochrone Object
******************************************

Isochrone 
=================
.. function:: class Isochrone(params)

   Create a theoretical cluster isochrone.

   :param logAge: The age of the isochrone, in log(years)
   :param AKs: The extinction of the isochrone, in magnitudes in the Ks filter
   :param distance: The distance to the isochrone, in parsecs
   :param metallicity: The metallicity of the isochrone, in [M/H]
   :param evo_models: The stellar evolution models used (:ref:`evo_models`)
   :param atm_func: The stellar atmosphere models used (:ref:`atmo_models`)
   :param red_law: The extinction law to be used (:ref:`ext_law`)
   :param mass_sampling: Sample the raw isochrone at every N
			 points, where N=mass_sampling. Must be an
			 integer value.
   :param wave_range: A list of length=2 that defines the min/max
		      wavelength of the stellar spectra.
   :param min_mass: Set the minimum mass of the isochone, in M_sun
   :param max_mass: Set the maximum mass of the isochone, in M_sun
   :param rebin: If true, rebin the resolution of the atmospheres so
		 they are the same as the Castelli+04 atmospheres


Example
---------


IsochronePhot
=================




