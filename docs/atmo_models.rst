.. _atmo_models:

========================================
Atmosphere Model Object
========================================
Stellar atmosphere models are defined as functions in
``popstar/atmospheres.py``. These can be called by::

 from popstar import atmospheres
 atmo = atmospheres.<function_name>

To call an atmosphere for a particular star, user must define the
metallicity ([Z]), temperature (in K), and gravity (in cgs)::

  spectrum = atmo(metallicity=0, temperature=5800, gravity=4)
  from astropy import units
  wave = spectrum.waveset.to(units.AA).value   # Wavelength in Angstroms
  flux = spectrum(spectrum.waveset)  # Flux (synphot Quantity; typically FLAM)

The atmosphere function is an input for the :ref:`isochrone_objects`,
and will automatically be used to define the
spectrum of each star in the isochrone model.

PopStar uses the stsynphot/synphot stack to extract the model atmosphere,
and the output spectrum is a synphot ``SourceSpectrum`` built from the CDBS grid via ``stsynphot.grid_to_spec``.

Below is a table of atmosphere model grids currently supported by
SPISEA. Note that the resolution column reports the original
resolution of the atmosphere model grid. These are available in the
`spisea_cdbs_highres.tar.gz` file on the installation page. However, the default grid SPISEA
uses has degraded the resolution of all atmosphere grids to R = 250
(the `spisea_cdbs.tar.gz` file). 

.. figure:: images/atm_models.png
	    :width: 900
            :height: 196
	    :align: center

Model Atmosphere Classes
-------------------------
.. autofunction:: atmospheres.get_merged_atmosphere

.. autofunction:: atmospheres.get_wd_atmosphere

.. autofunction:: atmospheres.get_bb_atmosphere
   
.. autofunction:: atmospheres.get_castelli_atmosphere

.. autofunction:: atmospheres.get_phoenixv16_atmosphere
		  
.. autofunction:: atmospheres.get_BTSettl_2015_atmosphere
		  
.. autofunction:: atmospheres.get_BTSettl_atmosphere
		 
.. autofunction:: atmospheres.get_wdKoester_atmosphere

.. autofunction:: atmospheres.get_kurucz_atmosphere

.. autofunction:: atmospheres.get_phoenix_atmosphere

