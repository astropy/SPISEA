.. _isochrone_objects:

================
Isochrone Object
================
The isochrone classes are defined in popstar/synthetic.py. The primary
inputs to a isochrone object are the stellar population age, distance,
total extinction, and metallicity, along with the :ref:`atmo_models`, 
:ref:`evo_models`, and :ref:`ext_law`. 

If the IsochronePhot sub-class is used then synthetic photometry
will be produced. The :ref:`filters` are defined as additional inputs.

An example of making an IsochronePhot object::

  from popstar import synthetic, evolution
  from popstar import atmospheres, reddening
  import numpy as np

  # Define isochrone input parameters
  logAge = np.log10(5*10**6.) # Age in log(years)
  dist = 4000 # distance in parsec
  AKs = 0.8 # extinction in Ks-band mags
  metallicity = 0 # Metallicity in [M/H]
  evo_model = evolution.MISTv1() 
  atm_func = atmospheres.get_merged_atmosphere
  red_law = reddening.RedLawHosek18b()

  # Specify filters for synthetic photometry. Here we 
  # use  the HST WFC3-IR F127M, F139M, and F153M filters
  filt_list = ['wfc3,ir,f127m', 'wfc3,ir,f139m', 
  'wfc3,ir,f153m']

  # Specify the directory we want the output isochrone
  # table saved in
  iso_dir = './isochrones/'

  # Make IsochronePhot object 
  my_iso = synthetic.IsochronePhot(logAge, AKs, dist, 
                            metallicity=0,
                            evo_model=evo_model, 
                            atm_func=atm_func,
                            red_law=red_law, 
                            filters=filt_list,
			    iso_dir=iso_dir)

See `Quick Start Example
<https://github.com/astropy/PyPopStar/blob/new_doc/docs/Quick_Start_Make_Cluster.ipynb>`_
for a detailed example for how to interact with the isochrone object.


Tips and Tricks: The IsochronePhot Object
-----------------------------------------

* It takes ~1-3 mins to make an IsochronePhot object for the first
  time. A FITS table is created with the stellar parameters and
  photometry for each star. This table is saved in the specified
  iso_dir, under the filename iso_<age>_<aks>_<dist>_<z>.fits.
  By default, iso_dir is set to the current working directory unless
  otherwise defined. 

  In future calls of IsochronePhot, the code will first check iso_dir
  to see if the appropriate FITS table already exists. If so,
  then it will simply read the table and be done. This saves
  significant amounts of computation time.

  * **WARNING**: When IsochronePhot checks to see if the desired
    isochrone table already exists, it checks all isochrone properties
    except for the photometric filters (evolution models, atmosphere
    models, and reddening law are encoded in the table meta-data).
    If any of these parameters do not match, then the isochrone will
    be re-calculated.

    However, to keep the isochrone filenames reasonable, only the
    age, extinction, distance, and metallicity are encoded in the
    filename itself. So, if the evolution model, atmosphere model, or
    reddening law have changed, the original file will be overwritten
    by the new isochrone.

    *To avoid files from being unintentially overwritten, we recommend
    that users specify different iso_dir paths when making isochrones
    with different evolution models, atmosphere models, or reddening
    laws.*
    
  * **WARNING**: IsochronePhot does not check existing
      isochrone tables to see if the photometric filters match
      those specified by the user. *So, if the user wishes to generate an
      isochrone with different filters, we recommend either using a
      different iso_dir path or setting the keyword recomp=True (see
      docs below).*

Base Isochrone Class
----------------------------
.. autoclass:: synthetic.Isochrone
	       :members: plot_HR_diagram, plot_mass_luminosity



Isochrone Sub-classes
-----------------------

.. autoclass:: synthetic.IsochronePhot
	       :show-inheritance:
		:members: make_photometry, plot_CMD, plot_mass_magnitude
