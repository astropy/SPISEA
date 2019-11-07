.. _getting_started:


==========================
Installation From Git
==========================
PyPopStar is hosted on GitHub `here
<https://github.com/astropy/PyPopStar>`_.
The master branch is the current release,
while the development branch is the development version.


.. _Dependencies:

Dependencies
----------------
* python (preferably via AstroConda, as it includes some necessary
  packages, like astropy and pysynphot)
* astropy
* pysynphot
* scipy
* numpy
* matplotlib
* STScI CDBS data package (Download from `here
  <http://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/synphot-throughput-tables.html>`_
  or directly via ftp `here
  <ftp://archive.stsci.edu/pub/hst/pysynphot>`_. You will need 2
  files: synphot1.tar.gz (63 MB unzipped) and synphot5.tar.gz (791 MB
  unzipped). See :ref:`set-up-cdbs` below for further instructions.)

.. _Build:

Build
------

Clone the git repo.
PyPopStar only consists of python codes, so you can either run the
standard ``python setup.py install``  or you can simply modify your path
environment variables to point to the git python location.

If the installation directory for PyPopStar (where you cloned the repo),
is ``path_to_pypopstar``, then for BASH, edit your ``.bash_profile`` to
include::

   export PYTHONPATH=$PYTHONPATH:<path_to_pypopstar>/PyPopStar

where you replace ``<path_to_pypostar>`` with the appropriate directory. 

For C-shell, edit your .cshrc to include::

    setenv PYTHONPATH <path_to_pypopstar>/PyPopStar

where you replace ``path_to_postar`` with the appropriate
directory. Note, the python path should point to the top-level "PyPopStar"
directory. The actual python functions should be found in the
``PyPopStar/popstar/`` sub-directory.

.. _set-up-cdbs:

Set Up CDBS Directory
-------------------
PyPopStar uses the STScI CDBS infrastructure to store
model atmospheres and HST filter functions. It expects
a top-level ``cdbs`` directory with 3 sub-directories: ``comp``, ``grid``,
and ``mtab``. The ``comp`` and ``mtab`` directories are obtained from the
``synphot1.tar.gz`` file you downloaded (grp/hst/cdbs/comp,
grp/hst/cdbs/mtab). The ``grid`` directory will eventually hold all of
the model atmospheres, but start by putting the Phoenix model
directory from ``synphot5.tar.gz`` there (grp/hst/cdbs/grid/phoenix).
You will add additional atmosphere model directories to ``cdbs/grid`` in
the :ref:`models` section below.

So, you should end up with a single ``cdbs`` directory with
the ``cdbs/comp``, ``cdbs/grid``, and ``cdbs/mtab`` directories
defined.

.. _models:

Stellar Evolution and Atmosphere Models
--------------------------------------
Once you have cloned the PyPopStar repository, you will need to download the
grids of evolution and atmosphere models. The evolution models are
stand-alone and specific to PyPopStar. The atmosphere models use the
STScI CDBS conventions and should be copied into your ``cdbs/grid`` directory.

The two files you need to download are:

* `popstar_models.tar.gz
  <http://astro.berkeley.edu/~jlu/popstar/popstar_models.tar.gz>`_. (2.6 GB; 18 GB unzipped)

* `postar_cdbs.tar.gz <http://astro.berkeley.edu/~jlu/popstar/popstar_cdbs.tar.gz>`_  (142 MB; 248 MB unzipped)

You can also optionally download a third file, which contains high-resolution versions of the atmospheres in ``popstar_cdbs.tar.gz``:

* `popstar_cdbs_highres.tar.gz <http://astro.berkeley.edu/~jlu/popstar/popstar_cdbs_highres.tar.gz>`_ (50 GB; 74 GB unzipped)

However, `popstar_cdbs_highres.tar.gz`` is not needed for PyPopStar
to run. PyPopStar uses the low-resolution atmospheres in
``popstar_cdbs.tar.gz`` for synthetic photometry by default, as
this is much faster and is sufficient in most applications. 

Once downloaded, ``popstar_cdbs.tar.gz`` (and
``popstar_cdbs_highres.tar.gz``, if desired) should be
expanded in  the directory that houses ``cdbs``.
The ``popstar_models.tar.gz`` file can be expanded
anywhere, as you will set an environment variable to its location.
However, we recommend downloading it to the same location as your cdbs directory. 

For instance, if your cdbs directory is ``/g/lu/models/cdbs/``
then you should put the .gz files in ``/g/lu/models/``
and then unzip them from there::

   cd /g/lu/models/
   tar xvf popstar_cdbs.tar.gz
   tar xvf popstar_models.tar.gz


``popstar_cdbs.tar.gz`` will put the model atmospheres in
``cdbs/grid``, and ``popstar_models.tar.gz`` will put the evolution
models in a directory called ``evolution/``

See
[README_EvolutionModels.txt](https://github.com/astropy/PyPopStar/blob/master/docs/README_EvolutionModels.txt)
for a description of the stellar evolution models and
[README_AtmoModels.txt](https://github.com/astropy/PyPopStar/blob/master/docs/README_AtmoModels.txt)
for a description of the stellar atmosphere models included in this download.
You can add additional evolution models and atmosphere models to PyPopStar if you wish.

.. _setup-paths:

Setup Path to Models
------------------

You need to notify python where the evolution and atmosphere models
live. This is done with two environment variables, ``PYSYN_CDBS`` and
``POPSTAR_MODELS``.  In your .cshrc or .bashrc, set them as the following::
  setenv PYSYN_CDBS /<path_to_models_directory>/cdbs
  setenv POPSTAR_MODELS /<path_models_directory>/


or::
  export PYSYN_CDBS=/<path_to_models_directory>/cdbs
  export POPSTAR_MODELS=/<path_to_models_directory>


Where the ``<path_to_models_directory>`` is where you downloaded and
expanded the ``popstar_cdbs.tar.gz`` and ``popstar_models.tar.gz``
files (e.g., where the ``cdbs/`` and ``evolution/`` directories live). 

.. _test-setup:

Testing PyPopStar Setup
---------------------------------------

If all goes well, you should be able to import any of the PyPopStar functions
in a python environment window using an import statement like those at the top
of the Quick Start Guide, such as::
    
    from popstar import synthetic
    
If the ``PYTHONPATH`` is broken, then python won't know where to find the popstar codes and
you will get an error reporting "No module named popstar". If the ``POPSTAR_MODELS`` or 
``PYSYN_CDBS`` paths are broken, then pypopstar won't know here to go to get the 
stellar models. and you will get either or both of the following warnings::

    UserWarning: PYSYN_CDBS is undefined; functionality will be SEVERELY crippled.
    
    UserWarning: POPSTAR_MODELS is undefined; functionality will be SEVERELY crippled.
    
Otherwise, you should be all set! You may get warnings that say::

    UserWarning: Extinction files not found in /u/mwhosek/models/cdbs/extinction
    
    UserWarning: No thermal tables found, no thermal calculations can be performed
    
But these can be safely ignored since PyPopstar doesn't use this functionality.

To further test PyPopstar, try running the :ref:`quick_start`. 
To test the range of evolution models, atmosphere models, and photometric
filters available in PyPopStar, run the ``test_evolution_models()``, ``test_atmospheres_models()``, and ``test_filters()`` functions in ``popstar/tests/test_models.py``. 
