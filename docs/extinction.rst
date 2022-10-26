.. _ext_law:

===============================
Extinction Law Object
===============================
The extinction law can be defined using the classes in ``spisea/reddening.py``. These can be called by::

  from spisea import reddening
  red_law = reddening.<redlaw_name>()

SPISEA uses the pysynphot framework to define the extinction law.
The output is a `pysynphot.reddening.CustomRedLaw
<https://pysynphot.readthedocs.io/en/latest/ref_api.html#module-pysynphot.extinction>`_
oject.
The reddening law is reported in terms of A_lambda / A_Ks, and thus is normalized to A_Ks = 1.

The red_law object is passed into the :ref:`isochrone_objects` in order to
define the extinction for the stars. See the `Quick Start <https://github.com/astropy/SPISEA/blob/main/docs/Quick_Start_Make_Cluster.ipynb>`_
for an example.

Available extinction laws:

* RedLawPowerLaw
* RedLawBrokenPowerLaw
* RedLawRiekeLebofsky
* RedLawCardelli
* RedLawIndebetouw05
* RedLawRomanZuniga07
* RedLawFitzpatrick09
* RedLawNishiyama09 (default)
* RedLawSchoedel10
* RedLawFritz11
* RedLawDamineli16
* RedLawDeMarchi16
* RedLawSchlafly16
* RedLawHosek18b
* RedLawNoguerasLara18
* RedLawNoguerasLara20


Extinction Law Classes
--------------------------

.. autoclass:: reddening.RedLawPowerLaw
	       :members: powerlaw

.. autoclass:: reddening.RedLawBrokenPowerLaw
	       :members: broken_powerlaw

.. autoclass:: reddening.RedLawRiekeLebofsky
	       :members: RiekeLebofsky85
			 
.. autoclass:: reddening.RedLawCardelli
	       :members: Cardelli89

.. autoclass:: reddening.RedLawIndebetouw05
	       :members: Indebetouw05

.. autoclass:: reddening.RedLawRomanZuniga07
	       :members: RomanZuniga07

.. autoclass:: reddening.RedLawFitzpatrick09
	       :members: Fitzpatrick09			 

.. autoclass:: reddening.RedLawNishiyama09
	       :members: Nishiyama09

.. autoclass:: reddening.RedLawSchoedel10
	       :members: Schoedel10

.. autoclass:: reddening.RedLawFritz11
	       :members: Fritz11
			 
.. autoclass:: reddening.RedLawDamineli16
	       :members: Damineli16

.. autoclass:: reddening.RedLawDeMarchi16
	       :members: DeMarchi16
			 
.. autoclass:: reddening.RedLawSchlafly16
	       :members: Schlafly16

.. _Hosek_new:
.. autoclass:: reddening.RedLawHosek18b
	       :members: Hosek18b

.. autoclass:: reddening.RedLawNoguerasLara18
	       :members: NoguerasLara18

.. autoclass:: reddening.RedLawNoguerasLara20
	       :members: NoguerasLara20


