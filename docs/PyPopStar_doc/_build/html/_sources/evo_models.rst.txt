.. _evo_models:

========================================
Evolution Model Object
========================================

Stellar evolution models are defined as classes in
popstar/evolution.py. These can be called by::

  from popstar import evolution
  evo = evolution.<model_name>()

The stellar evolution model classes currently supported by PyPopStar
are:

* MergedBaraffePisaEkstromParsec (default)
* MISTv1
* MergedPisaEkstromParsec
* MergedSiessGenevaPadova
* Baraffe15
* Ekstrom12
* Parsec
* Pisa

<Include table with model properties here>
  
The properties of the individual classes are detailed below.

.. _merged:

MergedBaraffePisaEkstromParsec
------------------------------------------------------
This is a combination of several different evolution models, and is
the default set used unless otherwise specified.

Component Models:

* Baraffe (`Baraffe et al. 2015 <https://ui.adsabs.harvard.edu//#abs/2015A&A...577A..42B/abstract>`_)
* Pisa (`Tognelli et al. 2011 <https://ui.adsabs.harvard.edu//#abs/2011A&A...533A.109T/abstract>`_)
* Geneva (`Ekstrom et al. 2012 <https://ui.adsabs.harvard.edu//#abs/2012A&A...537A.146E/abstract>`_)
* Parsec (version 1.2s, `Bressan+12 <https://ui.adsabs.harvard.edu//#abs/2012MNRAS.427..127B/abstract>`_)

For logAge < 7.4:

* Baraffe: 0.08 - 0.4 M_sun
* Baraffe/Pisa transition: 0.4 - 0.5 M_sun 
* Pisa: 0.5 M_sun to the highest mass in Pisa isochrone (typically 5 - 7 Msun)
* Geneva: Highest mass of Pisa models to 120 M_sun

For logAge > 7.4:
* Parsec v1.2s: full mass range

MISTv1
-------
Mesa Isochrone and Stellar Tracks (MIST) (`Choi et al. 2016 <https://ui.adsabs.harvard.edu//#abs/2016ApJ...823..102C/abstract>`_)

Versions 1.0, 1.2 supported. v1.2 is the default

Downloaded via web interpolator: http://waps.cfa.harvard.edu/MIST/interp_isos.html


MergedPisaEkstromParsec
-----------------------------------------
Same blend of models as :ref:`merged`, but without Baraffe+15.

MergedSiessGenevaPadova
-----------------------------------------
This is a combination of several different evolution models.

Component Models:

* Siess (`Siess et al. 2000 <https://ui.adsabs.harvard.edu//#abs/2000A&A...358..593S/abstract>`_)
* Geneva (`Meynet & Maeder 2003 <https://ui.adsabs.harvard.edu//#abs/2003A&A...404..975M/abstract>`_)
* Padova (`Marigo et al. 2008 <https://ui.adsabs.harvard.edu/abs/2008A%26A...482..883M/abstract>`_)

For logAge < 7.4:

* Siess: 0.1 - 7 M_sun
* Siess/Geneva transition: 7 - 9 M_sun
* Geneva: > 9 M_sun

For logAge > 7.4:

* Padova: full mass range

Baraffe15
--------------
Models published in `Baraffe et al. 2015 <https://ui.adsabs.harvard.edu//#abs/2015A&A...577A..42B/abstract>`_,
Downloaded from http://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/BHAC15_tracks

Esktrom12
---------------
Models published in `Ekstrom et al. 2012
<https://ui.adsabs.harvard.edu//#abs/2012A&A...537A.146E/abstract>`_,
Downloaded from http://obswww.unige.ch/Recherche/evoldb/index/Isochrone/

Parsec
------
Models published in `Bressan et al. 2012
<https://ui.adsabs.harvard.edu//#abs/2012MNRAS.427..127B/abstract>`_,
version 1.2s.
Downloaded from http://stev.oapd.inaf.it/cgi-bin/cmd

Parameters used in download:

* n_Reimers parameter (mass loss on RGB) = 0.2
* photometric system: HST/WFC3 IR channel
* bolometric corrections OBC from Girardi+08, based on ATLAS9 ODFNEW models
* Carbon star bolometric corrections from Aringer+09
* no dust
* no extinction
* Chabrier+01 mass function

Pisa
----
Models published in `Tognelli et al. 2011< https://ui.adsabs.harvard.edu//#abs/2011A&A...533A.109T/abstract>`_,
Downloaded from http://astro.df.unipi.it/stellar-models/index.php?m=1

Parameters used in download:

* Y = middle value of 3 provided (changes for different metallicities)
* mixing length = 1.68
* Deuterium fraction: 2*10^-5 for Z = 0.015, 0.03; 4*10^-4 for 0.005
