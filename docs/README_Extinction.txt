The extinction law can be defined using the classes in popstar/reddening.py. These can be called by:

> from popstar import reddening
> red_law = reddening.<redlaw_name>()

PopStar uses the pysynphot framework to define the extinction law. The output is a pysynphot.reddening.CustomRedLaw object (https://pysynphot.readthedocs.io/en/latest/ref_api.html#module-pysynphot.extinction). The reddening law is reported in terms of A_lambda / A_Ks, and thus is normalized to A_Ks = 1.

The red_law object is passed into the Isochrone object in order to define the extinction for all of the stars. See Quick_Start_Make_Cluster for an example.

Available extinction laws:

-RedLawPowerLaw
-RedLawCardelli
-RedLawDamineli16
-RedLawDeMarchi16
-RedLawFitzpatrick09
-RedLawFritz11
-RedLawHosek18
-RedLawHosek18b
-RedLawNishiyama09 (default)
-RedLawNoguerasLara
-RedLawRomanZuniga07
-RedLawRiekeLebofsky
-RedLawSchlafly16


============================
RedLawPowerLaw
============================
This define a power-law extinction law (A_lambda /propto lambda^-alpha).

Inputs:
alpha -- The power law exponent
K_wave -- Define the specific Ks wavelength to normalize the law too (in angstoms) 
wave_min, wave_max -- The shortward and longward bound of the extinction law (in angstroms)

Example:
> red_law = reddening.RedLawPowerLaw(2.21, 2.12, wave_min=0.8, wave_max = 3.0)

This defines the reddening law as A_lambda / A_Ks = lambda^-2.12, where the Ks wavelength is taken to be 2.12 microns. The law is defined from 0.8 -- 3.0 microns.

============================
RedLawCardelli
============================
Defines the Cardelli et al. 1989 law (https://ui.adsabs.harvard.edu//#abs/1989ApJ...345..245C/abstract)

Wave range: 0.5 - 3.0 microns

Inputs:
Rv -- The Cardelli+89 R_v parameter (R_v = 3.1 for average Milky Way ISM)


============================
RedLawDamineli16
============================
Defines the Damineli et al. 2016 law measured for Wd1 (https://ui.adsabs.harvard.edu//#abs/2016MNRAS.463.2653D/abstract)

Wave range: 0.5 - 8.0 microns

============================
RedLawDeMarchi16
============================
Defines De Marchi et al. 2016 law measured for 30 Dorodus (https://ui.adsabs.harvard.edu//#abs/2016MNRAS.455.4373D/abstract)

Wave range: 0.3 - 8 microns

============================
RedLawFitzpatrick09
============================
Defines Fitzpactrick et al. 2009 law (https://ui.adsabs.harvard.edu//#abs/2009ApJ...699.1209F/abstract)

Wave range: 0.3 - 3 microns

Inputs:
alpha -- free parameter in law (see their Eqn 5)
RV -- free parameter in law (see their Eqn 5)

============================
RedLawFritz11
============================
Defines Fritz et al. 2011 law for the Galactic Center (https://ui.adsabs.harvard.edu//#abs/2011ApJ...737...73F/abstract)

Wave range: 1.0 -- 19 microns

============================
RedLawHosek18
============================
Defines Hosek et al. 2018 law for the Arches cluster and Wd1 (https://ui.adsabs.harvard.edu//#abs/2018ApJ...855...13H/abstract)

Wave range: 0.7 - 3.54 microns

NOTE: this law has been updated (Hosek et al. 2019, Appendix B, https://ui.adsabs.harvard.edu//#abs/2019ApJ...870...44H/abstract), so should use the updated one instead (RedLawHosek18b, see below)

============================
RedLawHosek18b
============================
Defines Hosek et al. 2019 law for the Arches and Wd1 (https://ui.adsabs.harvard.edu//#abs/2019ApJ...870...44H/abstract). This should be used rather than RedLawHosek18.

Wave range: 0.7 - 3.54 microns

============================
RedLawNishiyama09
============================
Defines Nishiyama et al. 2009 law towards the center of the Milky Way (https://ui.adsabs.harvard.edu//#abs/2009ApJ...696.1407N/abstract). This is the default extinction law. 

Wave range: 0.3 -- 8.0 microns

============================
RedLawNoguerasLara18
============================
Defines Nogueras-Lara et al. 2018 law for the Galactic Center (https://ui.adsabs.harvard.edu//#abs/2018A&A...610A..83N/abstract)

Wave range: 0.8 - 2.5 microns

============================
RedLawRomanZuniga07
============================
Defines Roman-Zuniga et al. 2007 law for dense cloud core Barnard 59 (https://ui.adsabs.harvard.edu//#abs/2007ApJ...664..357R/abstract)

Wave range: 1.0 - 8.0 microns

============================
RedLawRiekeLebofsky
============================
Defines the Rieke & Lebofsky 1985 law for the Galactic Center (https://ui.adsabs.harvard.edu//#abs/1985ApJ...288..618R/abstract)

Wave range: 1 - 13 microns


============================
RedLawSchlafly16
============================
Defines the Schlafly et al. 2016 law (https://ui.adsabs.harvard.edu//#abs/2016ApJ...821...78S/abstract)

Wave range: 0.5 - 8 microns

Inputs:
AH_AKs -- Ratio of A_H / A_Ks, which sets the normalization of the law (see their Appendix)
x -- free parameter in extinction law (see their Eqn 6)
