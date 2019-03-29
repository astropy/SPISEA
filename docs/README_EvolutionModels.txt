============================
Authors
Matthew Hosek Jr, Casey Lam
============================
The stellar evolution model classes currently supported by PopStar
are:

-MergedBaraffePisaEkstromParsec (default)
-MISTv1
-MergedPisaEkstromParsec
-MergedSiessGenevaPadova
-Baraffe15
-Ekstrom12
-Parsec
-Pisa

Properties of these classes are described below.

============================
MergedBaraffePisaEkstromParsec
============================
This is a combination of several different evolution models, and is
the default set used unless otherwise
specified.

mass range: 0.08 - 120 M_sun
log(Age) range: 6.0 - 10.09
metallicity: solar
rotation: yes (Ekstrom only)

Composition:
-Baraffe et al. 2015
(https://ui.adsabs.harvard.edu//#abs/2015A&A...577A..42B/abstract)
-Pisa (Tognelli et al. 2011, https://ui.adsabs.harvard.edu//#abs/2011A&A...533A.109T/abstract)
-Geneva (Ekstrom et al. 2012,
https://ui.adsabs.harvard.edu//#abs/2012A&A...537A.146E/abstract)
-Parsec (version 1.2s, Bressan+12:
https://ui.adsabs.harvard.edu//#abs/2012MNRAS.427..127B/abstract)

For logAge < 7.4:
Baraffe+15: 0.08 - 0.4 M_sun
Baraffe+15/Pisa transition: 0.4 - 0.5 M_sun 
Pisa: 0.5 M_sun to the highest mass in Pisa isochrone (typically
5 - 7 Msun)
Geneva: Highest mass of Pisa models to 120 M_sun

For logAge > 7.4:
Parsec v1.2s: full mass range

=====================================
MISTv1
=====================================
-Mesa Isochrone and Stellar Tracks (MIST), Choi et al. 2016
(https://ui.adsabs.harvard.edu//#abs/2016ApJ...823..102C/abstract)
-versions 1.0, 1.2 supported. v1.2 is the default
-Downloaded via web interpolator: http://waps.cfa.harvard.edu/MIST/interp_isos.html

mass range: 0.1 - 300 M_sun
log(Age) range: 6.0 - 10.01
metallicity: solar (Apslund+09)
rotation: yes

============================
MergedPisaEkstromParsec
============================
Same blend of models as MergedBaraffePisaEkstromParsec, but without Baraffe+15.

mass range: 0.2 - 300 M_sun
log(Age) range: 6.0 - 8.0
metallicity: solar 
rotation: yes (Ekstrom only)

============================
MergedSiessGenevaPadova
============================
This is a combination of several different evolution models.

mass range: 0.1 - 150 M_sun
log(Age) range: 5.5 - 9.78
metallicity: solar 
rotation: no

Composition:
-Siess et al. 2000
(https://ui.adsabs.harvard.edu//#abs/2000A&A...358..593S/abstract)
-Geneva (Meynet & Maeder 2003,
https://ui.adsabs.harvard.edu//#abs/2003A&A...404..975M/abstract)
-Padova (CITE??)


For logAge < 7.4:
Siess+00: 0.1 - 7 M_sun
Interpolation between Siess+00 and Meynet & Maeder 2003: 7 - 9 M_sun
Meynet & Maeder 2003: > 9 M_sun

For logAge > 7.4:
Padova over full range

=================================
Baraffe15
=================================
Reference: Baraffe et al. 2015 (https://ui.adsabs.harvard.edu//#abs/2015A&A...577A..42B/abstract)
-Downloaded from http://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/BHAC15_tracks

mass range: 0.07 - 1.4 M_sun
log(Age) range: 5.7 - 10.0
metallicity: solar (Apslund+09)
mixing length: 1.6*pressure scale height
rotation: no

=================================
Esktrom12
=================================
Reference: Ekstrom et al. 2012
(https://ui.adsabs.harvard.edu//#abs/2012A&A...537A.146E/abstract)
-Downloaded from
http://obswww.unige.ch/Recherche/evoldb/index/Isochrone/


mass range: 0.8 - 300 M_sun
log(Age) range: 6.0 - 8.0
metallicity: Z = 0.14
rotation: both yes (0.4*v_crit) and no are available. Default = yes

========================
Parsec
========================
Version 1.2s
Reference: Bressan et al. 2012
(https://ui.adsabs.harvard.edu//#abs/2012MNRAS.427..127B/abstract)
-Downloaded from http://stev.oapd.inaf.it/cgi-bin/cmd

Parameters used:
-n_Reimers parameter (mass loss on RGB) = 0.2
-photometric system: HST/WFC3 IR channel
-bolometric corrections OBC from Girardi+08, based on ATLAS9 ODFNEW models
-Carbon star bolometric corrections from Aringer+09
-no dust
-no extinction
-Chabrier+01 mass function

mass range: 0.1 - 65 M_sun
log(Age) range: 6.6 - 10.12
metallicity: Z = 0.015, 0.04, 0.0005 (only solar supported thus far)
alpha enhancement: no
rotation: no

============================
Pisa
============================
Reference: Tognelli et al. 2011 https://ui.adsabs.harvard.edu//#abs/2011A&A...533A.109T/abstract
Isochrones downloaded from website:
http://astro.df.unipi.it/stellar-models/index.php?m=1

Parameters:
-Y = middle value of 3 provided (changes for different metallicities)
-mixing length = 1.68
-Deuterium fraction: 2*10^-5 for Z = 0.015, 0.03; 4*10^-4 for 0.005


mass range: 0.2 - 7 M_sun (upper mass decreases with age)
log(Age) range: 6.0 - 8.0
metallicity: Z = 0.015, 0.03, 0.005 (only solar supported thus far)
rotation: no
