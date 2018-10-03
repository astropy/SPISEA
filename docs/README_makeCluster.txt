#################################################################################
This file describes how to generate a cluster using the popstar package. 

Documentation Authors:
Matt Hosek
Siyao Jia

Notes/to-do bracketed by ======
#################################################################################

To generate a cluster, the user defines a cluster object which takes an isochone 
object, an IMF object, and a total cluster mass. Within the isochrone object, the 
user defines the cluster's age, distance, extinction, and metallity, as well as the 
evolution models, atmosphere models, and reddening law used. Within the IMF 
object, the user defines the functional form of the IMF used and the multiplicity 
parameters (default = no multiplicity). 

All code paths provided below have base path Popstar/popstar/

Step-by-step Instructions:

User can choose own models/reddening law to use, or they could just use 
the defaults:

1) Choose the stellar evolution model object that you would like to use.
code: evolution.py

example: evo_model = evolution.MergedBaraffePisaEkstromParsec

Options: 

a) MergedBaraffePisaEkstromParsec (default)
age range: logAge = 6.00 -- 9.99
mass range: 0.08 M_sun -- 500 M_sun (though age dependent)
metallicity: solar

For logAge < 7.4:
Baraffe et al. 2015 models from 0.08 - 0.4 M_sun
Transition region between 0.4 - 0.5 M_sun between Baraffe+15 and Pisa models
Pisa models (Tognelli et al. 2011) from 0.5 M_sun to the highest mass
available at that age (typically 5-7 M_sun), and then switching to
Geneva models (Ekstrom et al. 2012) for the rest of the way.

For logAge > 7.4:
Use the Parsec v1.2s (Bressan et al. 2012) models throughout entire
mass range

b) MergedPisaEkstromParsec: 
age range: logAge = 6.00 -- 8.0
mass range: 0.08 M_sun -- 500 M_sun (though age dependent)

Same as MergedBaraffePisaEkstromParsec, but without the Baraffe+15
models. As a result, the lower mass limit for logAge < 7.4 models is
0.2 M_sun


c) 


1) Choose the stellar atmosphere model object 
code: atmospheres.py

example: atmo_model = atmospheres.get_merged_atmosphere

Options:

a) get_merged_atmosphere (default)
Teff range: 1200 K -- ~50,000 K

BTStettl [1200 K -- 3200 K
PHOENIX v16 (Husser+2013) [2300 K - 5000 K]
PHOENIX/ATLAS merge [5000 K - 5500 K]
ATLAS (Castelli+2004) [5500K - 50000 K]



-atmosphere.get_merged_atmosphere
-reddening.RedLawNishiyama09






==============================
Matt 12/10/15
Functional running on own?
=============================
c) Geneva: Geneva models w/ rotation (Meynet + Maeder 2003)

d) Ekstrom: Geneva models w/ and w/o rotation (Ekstrom+2012)

e) Parsec: Parsec v1.2s models (Bressan+2012, Chen+2014)

f) Pisa: Pisa  models (Tognelli+2011)

g) Baraffe: Baraffe+2015 models (Baraffe+2015)

To Do: Add merged class for Lu+2013 models (Geneva + Seiss)?
================================

Atmosphere Models:
-code: atmospheres.py

Options:
a) get_merged_atmosphere: merge between 


b) get_kurucz_atmosphere: Kurucz+1993

c) get_castelli_atmosphere: Castelli+2004

d) get_nextgen_atmosphere: 

e) get_amesdusty_atmosphere: Allard+2000

f) get_phoenix_atmosphere: PHOENIX BT-Settl, Allard+2011

g) get_cmfgenRot_atmosphere, get_cmfgenNoRot_atmosphere: CMFGEN grid from Fierro+2015

h) get_phoenixv16_atmosphere: PHOENIX v16 grid from Husser+2013

i) get_atlas_phoenix_atmosphere: merge between ATLAS (Castelli+2004) and PHOENIX (Husser+2013)
between 5000 K - 5500 K

Reddening Law:
-code: reddening.py

Options:
a) RedLawNishiyama09: Nishiyama+2009, extinction law for Galactic Center [V - 13 microns]

b) RedLawCardelli: Cardelli+1989 (Rv = 3.1 default), extinction law for local ISM [V - 8 microns]

c) RedLawRomanZuniga07: RomanZuniga+2007, extinction law in Barnard 59 [1 - 13 microns]

d) RedLawRiekeLebofsky: Rieke + Lebofsky 1985, extinction law for Galactic Center [1 - 13 microns]


1) Create the cluster isochrone object.
-code: synthetic.py --> Isochrone, IsochronePhot

Options:

a) Isochone: returns the stellar parameters + spectra of the stars in the desired isochrone. 
DOES NOT SAVE OUTPUT!
Input: 
       age [log years]
       extinction [A_Ks mags]
       distance [pc]
       mass sampling of isochrone (default = 1)
       stellar evolution object (default = see step 0)
       stellar atmosphere function (default = see step 0)
       reddening law object (default = see step 0)
       wavelength range of output spectra (default = 5000 - 42500 Angstroms)

b) IsochronePhot: returns the stellar parameters + photometry in user-defined filters for the 
stars in the desired isochrone. SAVES OUPUT IN ISOCHRONE DIRECTORY DEFINED BY USER!
Input: 
       age [log years]
       extinction [A_Ks mags]
       distance [pc]
       isochrone directory (default = './')
       mass sampling of isochrone (default = 1)
       stellar evolution object	  (default = see step 0)
       stellar atmosphere function (default = see step 0)
       reddening law object (default = see step 0)
       filters used

Filters supported:
HST: all ACS and WFC3 filters supported 
(must use appropriate pysynphot syntax)
e.g. 'wfc3,ir,f125w'

NIRC2: 
'nirc2,J'
'nirc2,H'
'nirc2,K'
'nirc2,Kp'
'nirc2,Ks'
'nirc2,Lp'
'nirc2,Hcont'
'nirc2,Kcont'
'nirc2,FeII'
'nirc2,Brgamma'

2MASS: 
'2mass,J'
'2mass,H'
'2mass,Ks'

VISTA:
'vista,Z'
'vista,Y'
'vista,J'
'vista,H'
'vista,Ks'

DECam: u, g, r, i, z, y
'decam,u'
'decam,g'
'decam,r'
'decam,i'
'decam,z'
'decam,y'

PS1: g, r, i, z, y, w
'ps1,g'
'ps1,r'
'ps1,i'
'ps1,z'
'ps1,y'
'ps1,w'

JWST: F090W, F164N, F212N, F323N, F466N
'jwst,F090W'
'jwst,F164N'
'jwst,F212N'
'jwst,F323N'
'jwst,F466N'

filterList = {filtstring1, filtstring2, etc}
	   -use filtstrings as defined above

filter code: filters.py
filter file location: $POPSTAR_MODELS/filters

2) Make multiplicity object, if desired.
-code: imf/multiplicity.py

Options:
a) Multiplicity_Unresolved: Defines properties of stellar companions. 
Input:
	Multiplicity Fraction: A * m^gamma, m is primary mass (default: 0.44 * m^0.51) 
	Companion Star Fraction: B * m^beta, m is primary mass (default: 0.5 * m^0.45)
	q: Mass ratio between binaries; Q^power (default: 0.01 < Q < 1, power = -0.4)


3) Make IMF object. Will pass multiplicity object into IMF object.
-code: imf/imf.py 

Options:
a) IMF_broken_powerlaw

b) IMFSalpeter1955: Salpeter 1955

c) Miller_Scalo_1979: Miller+Scalo 1979

d) Kennicutt_1983: Kennicutt 1983

e) Kroupa_2001: Kroupa 2001

f) Weidner_Kroupa_2004: Weidner+Kroupa 2004

4) Make Cluster object
-code: synthetic.py

Options:
a) ResolvedCluster: Returns cluster with resolved stars or stellar systems (if binaries turned on). Parameters
and observables are returned for the individual stars.

b) UnresolvedCluster: Integrated spectrum of cluster over stellar population. 

