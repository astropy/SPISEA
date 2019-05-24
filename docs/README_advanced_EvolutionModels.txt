Advanced documentation for evolution models.

==========================
Baraffe15
==========================
How to make the isochrones from downloaded tracks:

0) SED voodoo magic
	-Need to reformat the downloaded tracks.dat fils so they
	can be read by python (specifically, "!" must be replaced by "#"). This is done using the following SED command:
		> sed 's/\!/#/' <oldfile> > <newfile>

1) popstar/evolution.py. class Baraffe15 --> tracks_to_isochrones
   	- Generate isochrones from evolutionary tracks. Resamples at 
	desired ages (6.0 < logAge 8.0, steps of 0.01)

	-additional function to test age-interpolated isochrone:
	class Baraffe15 --> test_age_interp

***File Locations:
-Original evolutionary tracks downloaded to
/g/lu/models/evolution/Baraffe15

-Isochrones stored in /g/lu/models/evolution/Baraffe15/iso

=============================
Ekstrom12
=============================
How to make the isochrones from downloaded tracks:

Note: website will let you download whatever time steps you want, but you are limited in the number of isochrones you can get per request. So, I had to download the 	isochrones in chunks and piece them together later

0) SED voodoo magic
	-Need to reformat the downloaded Isochrone*.dat files so they can be read by python (specifically, header lines need to be commented out). This is done using the 	following SED command:
		> sed 's/Isochrone/#Isochrone/' <oldfile> > <newfile>

1) popstar/evolution.py --> class EkstromStellarEvolution, function = create_iso
	-takes the Isochrones*.dat files downloaded online and puts them into a master iso.dat file with the proper format
		-will create tmp_*.dat files, one for each Isochrones*.dat file, which are easy to combine by hand in aquamacs

2) popstar/evolution.py --> class EkstromStellarEvolution, function = format_isochrones
	-parse the iso.dat file from create_iso function, create individual isochrone file for each given age
	
***File Locations

-Original downloaded files from website, plus tmp*.dat files after SED processed, plus master iso.dat file: /g/lu/models/evolution/Ekstrom2012/iso/z014

-rotating models: /g/lu/models/evolution/Ekstrom2012/iso/z014/rot
-non-rotating models: /g/lu/models/evolution/Ekstrom2012/iso/z014/norot

-evolutionary tracks in /g/lu/models/evolution/Ekstrom2012/ not used

=================================
Parsec
=================================
How to make the isochrones from downloaded tracks:

1) popstar/evolution.py --> class ParsecStellarEvolution, function format_isochrone
	-takes output*.dat file downloaded from web and organize into popstar format. Creates separate *.dat file for each age
	-can handle multiple metallicities at once, if desired

NOTE: no need to interpolate Parsec models ourselves since we can download whatever time resolution we want from the website

***File Locations:

-solar-metallicity isochrones: /g/lu/models/evolution/ParsecV1.2s/iso/z015
-sub-solar metallicity isochrones (not interpolated yet): /g/lu/models/evolution/ParsecV1.2s/iso/z005
-super-solar metallicity isochrones (not interpolated yet): /g/lu/models/evolution/ParsecV1.2s/iso/z04

-evolutionary tracks in /g/lu/models/evolution/ParsecV1.2s/tracks are not used


==================================
Pisa
==================================
How to make the isochrones from downloaded tracks:

1) popstar/evolution.py,  class PisaEvolutionaryModels, function = format_isochrones
	-rename isochrones downloaded from Pisa site to fit popstar iso*.dat scheme
	-can process multiple metallicities in same loop

2) popstar/evolution.py,  class PisaEvolutionaryModels, function = make_isochrone_grid
	-creates interpolated grid with the user-defined time sampling. Builds off of the online grid downloaded
		-hard-coded to create grid with delta-logAge = 0.01
		-only supports up to logAge = 8.0 (since this is the max age of the online grid)
		-can handle different metallicities if desired
	-this is a wrapper function which calls popstar/evolution.py ---> make_isochrone_pisa_interp
		-this in turn calls evolution.py --> get_orig_pisa_isochrones
	
	-paths are hard-coded, so don't need to be in any particular working directory here. Puts interpolated models in /g/lu/models/evolution/Pisa2011/iso/<metallicity>

NOTE: By design, this process interpolates over isochrones, rather than evolutionary tracks. Testing showed that interpolation over the isochrones was more accurate than interpolating over the isochrones

***File Locations:

-Original downloaded isochrones (needed for interpolation, DO NOT DELETE): /g/lu/models/evolution/Pisa2011/iso/
	-metallicity subdirectory: z015

-Interpolated solar metallicity isochrones: /g/lu/models/evolution/Pisa2011/iso/z015

-Non-interpolated non-solar metallicity isochrones: /g/lu/models/evolution/Pisa2011/iso/z005, z03

===================================
MISTv1
===================================
How to make the isochrones from downloaded tracks:

1) popstar/evolution.py --> class MISTv1, function format_isochrone
   -takes MIST_iso* files downloaded from web and organize into popstar format. Creates separate *.dat file for each age
   -can handle multiple metallicities at once, if desired

Note: No need to interpolate MIST isochrones since we can download them at any time step we want

***File Locations:
Original files are located in appropiate metallicity subdirectory under /g/lu/models/evolution/MISTv1

***Gory details on how to download isochrones below -CYL

1) Go to MIST isochrone page [http://waps.cfa.harvard.edu/MIST/interp_isos.html]
   Note: We have been downloading isochrones with the following convention:
         Version: 1.2
         Rotation: Initial v/v_crit = 0.4
         Age: Log10 Scale with spacing of 0.01
         Composition: Fe/H = 0
         Output Option: Theoretical, Full (79 columns)
    There's a limit on how many isochrones you can download at a time, so break up downloads into smaller chunks if necessary. 
    
2) Move the downloaded .iso file into the directory where you want to unpack the isochrone file.
   Note: For our current purposes this is /g/lu/models/evolution/MISTv1/v1.2/iso/z015. 
   But you should put it into the directory of the relevant version/metallicity.

3) Unpack the isochrone with format_isochrones function. An example of how to call in iPython:
   ```
   from popstar import evolution
   import numpy as np
   evo = evolution.MISTv1()
   input_iso_dir = /g/lu/models/evolution/MISTv1/v1.2/iso
   metallicity_list = np.array(['z015'])
   evo.format_isochrones(input_iso_dir, metallicity_list)
   ```
