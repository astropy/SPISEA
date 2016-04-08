====================
PopStar
====================

INSTALL (from git)

* python (preferably via ureka, as it includes some necessary packages, like astropy and pysynphot)
    python dependencies:
    
    * astropy
    * pysynphot
    * ATpy (https://pypi.python.org/pypi/ATpy) 
    * asciitable (http://www.stecf.org/software/PYTHONtools/astroasciidata/) - note: if python was installed via ureka, the asciitable install folder may need to be copied from the default python site-packages folder to ureka's site-packages folder
    
* Github respository dependencies (git clone these)
    * jluastro/JLU-python-code 
    * jluastro/nirc2

Add all of these to your python path in your .cshrc::

    set CODE_DIR = /path/to/top/level/python/code
    setenv PYTHONPATH ${CODE_DIR}/popstar/PopStar:${CODE_DIR}/jlu/JLU-python-code-master:${CODE_DIR}/jlu/JLU-python-code-master/jlu/gc:${CODE_DIR}/nirc2:${CODE_DIR}/ATpy-0.9.7:${CODE_DIR}/asciidata/asciidata-1.1.1

Once you have cloned the popstar repository, you will need to download the
grid of evolution and atmosphere models. The evolution models can be
found here::

    GIT SITE -- load these up
    /g/lu/models/evolution/

Ready:
 * * pisa_ekstrom_parsec
 * * pisa_ekstrom_parsec/norot

Data files ready but not yet supported in code:
 * Ekstrom2012 (rot and norot)
 * ParsecV1.2S
 * Pisa2011
 * geneva
 * merged
 * * siess_meynetMaeder_padova
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

In your .cshrc or .bashrc, set the `PSYN_CDBS` variable to point to the
CDBS directory:

.. highlight:: c

    setenv PSYN_CDBS /g/lu/models/cdbs
    export PSYN_CDBS=/g/lu/models/cdbs


Fubar2

.. _Astropy: http://www.astropy.org/
.. _git: http://git-scm.com/
.. _github: http://github.com
.. _Cython: http://cython.org/
