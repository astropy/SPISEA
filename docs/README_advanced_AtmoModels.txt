Advanced documentation for atmosphere models.

========================
Phoenixv16 Models
========================
Paper reference: Husser et al. 2013, A&A, 553, A6
downloaded via ftp from website: http://phoenix.astro.physik.uni-goettingen.de/?page_id=15

Copy wave file into each metallicity subdirectory. E.g.: `PHOENIX-ACES-AGSS-COND-2011_AtmosFITS_Z-0.0/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits`

### Code procedure:
1) popstar/atmospheres.py --> organize_PHOENIXv16_atmospheres
	- organizes files into subdirectory for the corresponding metallicity
	
	```py
	from popstar import atmospheres
	atmospheres.organize_PHOENIXv16_atmospheres('/u/abhimat/models/atmospheres/PHOENIXv16_Husser+13/PHOENIX-ACES-AGSS-COND-2011_AtmosFITS_Z-1.0/', met_str='m10')
	```
	
	- Copy metallicity subdirectories into cdbs grid directory for subsequent steps. E.g.: copy into '/u/abhimat/models/cdbs/grid/phoenix_v16/'

2) popstar/atmospheres.py --> make_PHOENIXv16_catalog
	- construct the cdbs catalog for the atmospheres organized into subdirectory for current metallicity
	- catalog.fits file is created in the current directory in which the script is run (e.g. '/u/abhimat/models/cdbs/grid/phoenix_v16/')
	
	```py
	from popstar import atmospheres
	
	atmospheres.make_PHOENIXv16_catalog('phoenixm10/', met_str='m10')
	```

3) popstar/atmospheres.py --> cdbs_PHOENIXv16
	- puts atmosphere files into cdbs format...important because we have a flux unit conversion here
	- needs to be run for every metallicity subdirectory
	
	```py
	from popstar import atmospheres
	atmospheres.cdbs_PHOENIXv16('phoenixm10/')
	```
	

4) popstar/atmospheres.py --> rebin_phoenixv16
	- Rebin phoenixv16 model atmospheres to same resolution as ck04 atlas, to enhance photometric performance
	```py
	from popstar import atmospheres
	atmospheres.rebin_phoenixV16('/u/abhimat/models/cdbs/')
	```

***File locations:
-raw models from online downloaded to /g/lu/model/atmospheres/PHEONIXv16_Husser+13
-raw models organized into phoenixm00 subdirectory (output from code step 1 above): /g/lu/model/atmospheres/PHEONIXv16_Husser+13
-cdbs-ready model atmospheres, high resolution: /g/lu/models/cdbs/grid/phoenix_v16
-cdbs-ready model atmospheres, low resolution: /g/lu/models/cdbs/grid/phoenix_v16_rebin

-old version of phoenix models: /g/lu/models/cdbs/grid/phoenix

===========================
CMFGEN Fierro+15 Grid
===========================
Paper Reference: 
Fierro et al. 2015, PASP, 127, 951

Models Downloaded from https://sites.google.com/site/fluxesandcontinuum/home
	-note: want the .flx models, which includes the spectral lines + continuum, not the .con models, which is only the continuum
	-two grids: one with and one without rotation (0.4 v_crit)
 
***Code procedure:
1) popstar/atmospheres.py --> download_CMFGEN_atmospheres
	-automated code to download bulk of models. NOTE: NOT 100% successful! Does miss a few models, user has to go back and add them by hand
	-note: requires tables describing the paramters of the rot and non-rot atmospheres. Paper doesn't provide electronic version of this, so I had to do this by hand. Decidedly 	not fun.

2) popstar/atmospheres.py --> organize_CMFGEN_atmospheres
	-organize atmospheres from downloaded state into cdbs-ready state. Creates two subdirectories: cmfgenF15_rot and cmfgenF15_norot

3) popstar/atmospheres.py --> make_CMFGEN_catalog
	-make cdbs catalog.fits for the cmfgen models. Must be run for rot and norot directories independently

4) popstar/atmospheres.py --> cdbs_cmfgen
	-take atmospheres in subdirectories created by organize_CMFGEN_atmospheres and convert them to cdbs format

5) popstar/atmospheres.py ---> rebin_cmfgen
	-rebin the cmfgen atmospheres into ck04 resolution, for improved photometric performance. "rot" flag in code is set to determine if code operates on rotating or non-	rotating codes

***File Locations:
-raw downloaded models (without flux normalization, which is what we want): /g/lu/models/atmospheres/CMFGEN_Fierro+15/ir/orig
	-with flux normalization can be found in CMFGEN_Fierro+15/flux_normalized_models/
-output from organize_CMFGEN_atmospheres (cmfgenF15_rot/norot): /g/lu/models/atmospheres/CMFGEN_Fierro+15/ir/
	-Table_rot/norot.txt files are in tese directories as well

-cdbs model atmospheres, original resolution: /g/lu/models/cdbs/grid/cmfgen_rot/norot
-cdbs model atmospheres, rebinned: /g/lu/models/cdbs/grid/cmfgen_rot/norot_rebin

===========================
BTSettl_CIFITS_2011_2015 Grid
===========================
Downloaded from https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/ on 9/28/16
-solar metallicity, no alpha enhancement

***Code procedure:
1) popstar/atmospheres.py --> organize_BTSettl_atmospheres 
    -organize atmospheres from downloaded state into cdbs-ready
    state. Moves them to /g/lu/model/atmospheres/BTSettl_2015

2) popstar/atmospheres.py --> make_BTSettl_catalog
	-construct the cdbs catalog for the atmospheres  

2.5) Need to make sure cdbs atmospheres have unique wavelength entries, otherwise pysynphot Icat load won't work. 
popstar/atmopsheres.py ---> make_wavelength_unique()

3) popstar/atmospheres.py --> rebin_BTSettl
	-rebin BTSettl model atmospheres to same resolution as ck04 atlas, to enhace photometric performance
	-these spectra moved to BTSettl_2015_rebin

***File Locations:
-raw models from online downloaded to /g/lu/model/atmospheres/BTSettl_CIFITS2011_2015
-cdbs-ready model atmospheres, high resolution: /g/lu/models/cdbs/grid/BTSettl_2015
-cdbs-ready model atmospheres, low resolution: /g/lu/models/cdbs/grid/BTSettl_2015_rebin
 





