.. _ext_law:

===============================
Extinction Law Object
===============================
The extinction law can be defined using the classes in popstar/reddening.py. These can be called by::

  from popstar import reddening
  red_law = reddening.<redlaw_name>()

PyPopStar uses the pysynphot framework to define the extinction law.
The output is a `pysynphot.reddening.CustomRedLaw
<https://pysynphot.readthedocs.io/en/latest/ref_api.html#module-pysynphot.extinction>`_
oject.
The reddening law is reported in terms of A_lambda / A_Ks, and thus is normalized to A_Ks = 1.

The red_law object is passed into the Isochrone object in order to define the extinction for all of the stars. See Quick_Start_Make_Cluster for an example.

Available extinction laws:

* RedLawPowerLaw
* RedLawCardelli
* RedLawDamineli16
* RedLawDeMarchi16
* RedLawFitzpatrick09
* RedLawFritz11
* RedLawHosek18 (deprecated)
* RedLawHosek18b
* RedLawNishiyama09 (default)
* RedLawNoguerasLara18
* RedLawRomanZuniga07
* RedLawRiekeLebofsky
* RedLawSchlafly16


RedLawPowerLaw
--------------------------
This define a power-law extinction law (A_lambda /propto lambda^-alpha).

Inputs:

alpha -- The power law exponent
K_wave -- Define the specific Ks wavelength to normalize the law too (in angstoms) 
wave_min, wave_max -- The shortward and longward bound of the extinction law (in angstroms)

Example::

  red_law = reddening.RedLawPowerLaw(2.21, 2.12, wave_min=0.8, wave_max = 3.0)

This defines the reddening law as A_lambda / A_Ks = lambda^-2.12, where the Ks wavelength is taken to be 2.12 microns. The law is defined from 0.8 -- 3.0 microns.


RedLawCardelli
-----------------------------
Defines the extinction law of `Cardelli et al. 1989 <https://ui.adsabs.harvard.edu//#abs/1989ApJ...345..245C/abstract>`_

Wave range: 0.5 - 3.0 microns

Inputs:

Rv -- The Cardelli+89 R_v parameter (R_v = 3.1 for average Milky Way ISM)


RedLawDamineli16
----------------------------
Defines the extinction law of `Damineli et al. 2016
<https://ui.adsabs.harvard.edu//#abs/2016MNRAS.463.2653D/abstract>`_
derived for Wd1

Wave range: 0.5 - 8.0 microns

RedLawDeMarchi16
----------------------------
Defines extinction law from `De Marchi et al. 2016
<https://ui.adsabs.harvard.edu//#abs/2016MNRAS.455.4373D/abstract>`_
derived for 30 Dorodus.

Wave range: 0.3 - 8 microns


RedLawFitzpatrick09
----------------------------
Defines the extinction law from `Fitzpactrick et al. 2009 <https://ui.adsabs.harvard.edu//#abs/2009ApJ...699.1209F/abstract>`_

Wave range: 0.3 - 3 microns

Inputs:

alpha -- free parameter in law (see their Eqn 5)
RV -- free parameter in law (see their Eqn 5)

RedLawFritz11
-----------------------------
Defines extinction law from `Fritz et al. 2011
<https://ui.adsabs.harvard.edu//#abs/2011ApJ...737...73F/abstract>`_
for the Galactic Center.

Wave range: 1.0 -- 19 microns

.. _Hosek_old:
RedLawHosek18
-----------------------------
Defines extinction law from `Hosek et al. 2018
<https://ui.adsabs.harvard.edu//#abs/2018ApJ...855...13H/abstract>`_
for the Arches Cluster and Wd1.

Wave range: 0.7 - 3.54 microns

NOTE: this law has revised to :ref:`Hosek_new`, which should be used
instead.

.. _Hosek_new:
RedLawHosek18b
-----------------------------
Defines extinction law from `Hosek et al. 2019
<https://ui.adsabs.harvard.edu//#abs/2019ApJ...870...44H/abstract>`_
for the Arches cluster adn Wd1. This should be used rather than :ref:`Hosek_old`

Wave range: 0.7 - 3.54 microns


RedLawNishiyama09
------------------------------
Defines extinction law from `Nishiyama et al. 2009
<https://ui.adsabs.harvard.edu//#abs/2009ApJ...696.1407N/abstract>`_
toward the Galactic Center. This is the default extinction law. 

Wave range: 0.3 -- 8.0 microns


RedLawNoguerasLara18
-------------------------------
Defines extinction law from `Nogueras-Lara et al. 2018
<https://ui.adsabs.harvard.edu//#abs/2018A&A...610A..83N/abstract>`_
for the Galactic Center.

Wave range: 0.8 - 2.5 microns


RedLawRomanZuniga07
-------------------------------
Defines extinction law from `Roman-Zuniga et al. 2007
<https://ui.adsabs.harvard.edu//#abs/2007ApJ...664..357R/abstract>`_
for the dense cloud core Barnard 59.

Wave range: 1.0 - 8.0 microns


RedLawRiekeLebofsky
------------------------------
Defines the extinction law from `Rieke & Lebofsky 1985
<https://ui.adsabs.harvard.edu//#abs/1985ApJ...288..618R/abstract>`_
for the Galactic Center

Wave range: 1 - 13 microns


RedLawSchlafly16
-----------------------------
Defines the extinction law from `Schlafly et al. 2016 <https://ui.adsabs.harvard.edu//#abs/2016ApJ...821...78S/abstract>`_

Wave range: 0.5 - 8 microns

Inputs:

AH_AKs -- Ratio of A_H / A_Ks, which sets the normalization of the law (see their Appendix)
x -- free parameter in extinction law (see their Eqn 6)
