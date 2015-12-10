######
This file describes how to generate a cluster using the popstar package. 

12/9/15
Matt Hosek
Lucy Jia
####

To generate a cluster, the user defines a cluster object which takes an isochone object, an IMF object, and a total cluster mass. Within the isochrone object, the user defines the cluster's age, distance, extinction, and metallity. Within the IMF object, the user defines the functional form of the IMF used and the multiplicity parameters (default = no multiplicity). 

Step-by-step Instructions:

1) Create the cluster isochrone object.
-code: synthetic.py --> Isochrone, IsochronePhot

-Isochrone object returns the stellar parameters + spectra of the stars in the desired isochrone. DOES NOT SAVE OUTPUT!
Input: 
       age (log years)
       extinction (A_Ks mags)
       distance (pc)
       mass sampling of isochrone
       stellar evolution object
       stellar atmosphere object
       reddening law object
       wavelength range of output spectra

-IsochronePhot object returns the stellar parameters + photometry in user-defined filters for the stars in the desired isochrone. DOES SAVE OUPUT!
Input: 
       age (log years)
       extinction (A_Ks mags)
       distance (pc)
       mass sampling of isochrone
       isochrone directory
       stellar evolution object
       stellar atmosphere object
       reddening law object
       filters used

2) Make multiplicity object, if desired.


3) Make IMF object. Will pass multiplicity object into IMF object.


4) Make Cluster object




