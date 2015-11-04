====================
PopStar
====================

INSTALL (from git)

Dependencies:
astropy
pysynphot

Once you have cloned the repository, you will need to download the
grid of evolution and atmosphere models. The evolution models can be
found here:

GIT SITE -- load these up
/g/lu/models/evolution/

Ready:
* Ekstrom2012 (rot and norot)
* ParsecV1.2S
* Pisa2011
* geneva
* merged
** pisa_ekstrom_parsec
** pisa_ekstrom_parsec/norot
** siess_meynetMaeder_padova
* padova
* pallaStahler1999
* siess2000

Update:
* STERN2011


The atmosphere models use the CDBS framework (which is supported by
pysynphot). We have added new grids of atmosphere models (including
merged sets of atmospheres). You can find the complete CDBS grid you
can download from here:

GIT SITE -- load these up
/g/lu/models/cdbs

Ready: 
* phoenix_v16 - high resolution (reformatted for CDBS)
* phoenix_v16_rebin - downgraded to improve synthetic photometry
performance.



MODEL LOCATION

You need to notify python where these models are going to live. This
is done in two steps.

In your .cshrc or .bashrc, set the PSYN_CDBS variable to point to the
CDBS directory:

.. highlight:: c

    setenv PSYN_CDBS /g/lu/models/cdbs
    export PSYN_CDBS=/g/lu/models/cdbs

Fubar


.. _Astropy: http://www.astropy.org/
.. _git: http://git-scm.com/
.. _github: http://github.com
.. _Cython: http://cython.org/
