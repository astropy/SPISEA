#################################################################################
This file describes how to generate a cluster using the popstar package. 

Authors:
Matt Hosek
Lucy Jia

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

0) Choose the stellar evolution models, stellar atmosphere models, and reddening 
law to be used in the isochrone.

User can choose own models/reddening law to use, or they could just use 
the defaults:
-evolution.MergedBaraffePisaEkstromParsec
-atmosphere.get_merged_atmosphere
-reddening.RedLawNishiyama09

If user does not define own model/reddening law, the defaults will automatically 
be adopted in the isochrone.

Evolution Models:
-code: evolution.py

Options:
a) MergedBaraffePisaEkstromParsec:


b) Merged PisaEkstromParsec: 


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
PHOENIX v16 (Husser+2013) [2300 K - 5000 K]
PHOENIX/ATLAS merge [5000 K - 5500 K]
ATLAS (Castelli+2004) [5500K - 50000 K]

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
HST: all filters supported in pysynphot (must use appropriate
pysynphot syntax)

NIRC2: J, H, K, Kp, Ks, Lp, Hcont, Kcont, FeII, Brgamma 

VISTA: Z, Y, J, H, Ks

DECam: u, g, r, i, z, y

PS1: g, r, i, z, y, w

JWST: F090W, F164N, F212N, F323N, F466N

filter syntax (for non-HST filters):
filterList = {'<user_name>':'<instrument>,<filter>'}
	   -<user_name> is user specified. The output magnitudes will
	   be in a column names mag_<user_name>
	   -<instrument>: name of telescope/instrument for the desired
	   filter. All letters lower case (i.e. vista, decam, ps1, etc)
	   -<filter>: name of filter, with proper case (see list above)

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

