.. _getting_started:


==========================
Install From Git
==========================
If you are downloading the code from scratch, please follow the
instructions below. If you had already downloaded version 1 of the
code and are switching to version 2, please see :ref:`version`. 

SPISEA is hosted on `GitHub <https://github.com/astropy/SPISEA>`_.
To begin, clone the git repository.
The ``main`` branch contains the current release,
while the ``dev`` branch is for code development.

.. _Dependencies:

Dependencies
----------------
* python (>=3.7; preferably via AstroConda, as it includes some necessary
  packages, like astropy and pysynphot)
* astropy
* pysynphot
* scipy
* numpy (>= v1.17)
* matplotlib
* STScI CDBS data package (Download from the `website
  <http://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/synphot-throughput-tables.html>`_
  or directly via ftp `here
  <ftp://archive.stsci.edu/pub/hst/pysynphot>`_. You will need 2
  files: synphot1.tar.gz (63 MB unzipped) and synphot5.tar.gz (791 MB
  unzipped). See :ref:`set-up-cdbs` below for further instructions.)

.. _Build:

Build
------
Once the repository is cloned, you need to add SPISEA to your
``PYTHONPATH`` environment variable. For example, if you are using
BASH, then you need to add the following line to your ``bash_profile``::
  
   export PYTHONPATH=$PYTHONPATH:<path_to_SPISEA>

where ``<path_to_SPISEA>`` is the location of the top-level SPISEA
directory. 

.. _set-up-cdbs:

Set Up CDBS Directory
---------------------------------
SPISEA uses the STScI CDBS infrastructure to store
model atmospheres and HST filter functions. It requires
a top-level ``cdbs`` directory with 3 sub-directories: ``comp``, ``grid``,
and ``mtab``. The ``comp`` and ``mtab`` directories are downloaded as
part of the ``synphot1.tar.gz`` file in the STSci CDBS download (grp/hst/cdbs/comp,
grp/hst/cdbs/mtab). The ``grid`` directory will eventually hold all of
the model atmospheres, but start by putting the Phoenix model
directory from ``synphot5.tar.gz`` there (grp/hst/cdbs/grid/phoenix).
You will add additional atmosphere model directories to ``cdbs/grid`` in
the :ref:`models` section below.

So, you should end up with a single ``cdbs`` directory with
the ``cdbs/comp``, ``cdbs/grid``, and ``cdbs/mtab`` sub-directories.

.. _models:

Stellar Evolution and Atmosphere Models
-------------------------------------------------------
SPISEA requires grids of stellar evolution and atmosphere models in
order to function. The evolution models are
stand-alone and specific to SPISEA. The atmosphere model grids use the
STScI CDBS conventions and should be placed in the ``cdbs/grid`` directory.

You will need to download 2 files:

* `popstar_models.tar.gz
  <http://astro.berkeley.edu/~jlu/popstar/popstar_models.tar.gz>`_. (2.6 GB; 18 GB unzipped)

* `postar_cdbs.tar.gz <http://astro.berkeley.edu/~jlu/popstar/popstar_cdbs.tar.gz>`_  (142 MB; 248 MB unzipped)

You can also optionally download a third file, which contains higher-resolution versions of the atmospheres in ``spisea_cdbs.tar.gz``:

* `popstar_cdbs_highres.tar.gz <http://astro.berkeley.edu/~jlu/popstar/popstar_cdbs_highres.tar.gz>`_ (50 GB; 74 GB unzipped)

SPISEA uses the low-resolution atmospheres (R = 250) in
``spisea_cdbs.tar.gz`` by default, as
it is then much faster to calculate synthetic photometry and
is sufficient in most applications. However, the user can change
this default; see  :ref:`atmo_models` for
more details. 

Once downloaded, ``spisea_cdbs.tar.gz`` (and
``spisea_cdbs_highres.tar.gz``, if desired) should be
expanded in  the ``cdbs`` directory. The output directories
will automatically be placed in the correct locations. 
The ``spisea_models.tar.gz`` file can be expanded
anywhere, though we recommend expanding it in the same location as 
your ``cdbs`` directory. 

For example, if your cdbs directory is ``/g/lu/models/cdbs/``
then you should put the .gz files in ``/g/lu/models/``
and then unzip them from there::

   cd /g/lu/models/
   tar xvf spisea_cdbs.tar.gz
   tar xvf spisea_models.tar.gz


``spisea_cdbs.tar.gz`` will put the model atmospheres in
``cdbs/grid``, and ``spisea_models.tar.gz`` will put the evolution
models in a directory called ``evolution/``

.. _setup-paths:

Set Paths to Models
-------------------------------------

You need to notify python where the evolution and atmosphere models
live. This is done by setting two environment variables, ``PYSYN_CDBS`` and
``SPISEA_MODELS``, to point to the ``cdbs`` and ``models``
directories (i.e. the directory where the ``evolution/`` directory
lives), respectively. For example, in your .bash_profile::
  
  export PYSYN_CDBS=/<path_to_cdbs_directory>
  export SPISEA_MODELS=/<path_to_models_directory>


.. _test-setup:

Testing Your SPISEA Setup
---------------------------------------

If all goes well, you should be able to import any of the SPISEA
functions an import statement like those at the top
of the Quick Start Guide, such as::
    
    from spisea import synthetic

You may get warnings that Extinction or thermal files are missing,
such as::

    UserWarning: Extinction files not found in /u/mwhosek/models/cdbs/extinction
    
    UserWarning: No thermal tables found, no thermal calculations can be performed
    
However, these can be safely ignored since SPISEA doesn't use those functionalities.

To further test your SPISEA install, try running the `Quick Start
notebook
<https://github.com/astropy/SPISEA/blob/master/docs/Quick_Start_Make_Cluster.ipynb>`_.
It is also located in ``SPISEA/docs``.

To test the full range of
evolution models, atmosphere models, and photometric filters,
run the test functions by going into the ``SPISEA/spisea`` directory and running::

    py.test

This will trigger the test functions we have implemented. If all is
well, you shouldn't get any errors. Warnings are (generally) fine.

TroubleShooting
-----------------------
If SPISEA is not properly in your ``PYTHONPATH``, then when you try
to import the SPISEA functions you will get an error message
reporting ``No module named spisea``.

If the ``SPISEA_MODELS`` or ``PYSYN_CDBS`` paths are broken, then
SPISEA won't know where to get the stellar models.
When trying to import ``spisea/synthetic.py``, you will get
either or both of the following warnings::

    UserWarning: PYSYN_CDBS is undefined; functionality will be SEVERELY crippled.
    
    UserWarning: SPISEA_MODELS is undefined; functionality will be SEVERELY crippled.
    
