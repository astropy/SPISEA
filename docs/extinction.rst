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
* RedLawCardelli
* RedLawDamineli16
* RedLawDeMarchi16
* RedLawFitzpatrick09
* RedLawHosek18 (deprecated)
* RedLawHosek18b
* RedLawNishiyama09 (default)
* RedLawNoguerasLara18
* RedLawRomanZuniga07
* RedLawRiekeLebofsky
* RedLawSchlafly16


  
Extinction Law Classes
--------------------------

.. autoclass:: reddening.RedLawPowerLaw
	       :members: powerlaw

.. autoclass:: reddening.RedLawBrokenPowerLaw
	       :members: broken_powerlaw
			 
.. autoclass:: reddening.RedLawCardelli
	       :members: Cardelli89

.. autoclass:: reddening.RedLawDamineli16
	       :members: Damineli16

.. autoclass:: reddening.RedLawDeMarchi16
	       :members: DeMarchi16
			 
.. autoclass:: reddening.RedLawFitzpatrick09
	       :members: Fitzpatrick09

.. autoclass:: reddening.RedLawFritz11
	       :members: Fritz11

.. _Hosek_old:
.. autoclass:: reddening.RedLawHosek18
	       :members: Hosek18

.. _Hosek_new:
.. autoclass:: reddening.RedLawHosek18b
	       :members: Hosek18b

.. autoclass:: reddening.RedLawNishiyama09
	       :members: Nishiyama09

.. autoclass:: reddening.RedLawNoguerasLara18
	       :members: NoguerasLara18

.. autoclass:: reddening.RedLawRomanZuniga07
	       :members: RomanZuniga07

.. autoclass:: reddening.RedLawRiekeLebofsky
	       :members: RiekeLebofsky85

.. autoclass:: reddening.RedLawSchlafly16
	       :members: Schlafly16
