.. _filters:

========================
Photometric Filters
========================

The user can specify what filters are used for synthetic photometry
when defining the :ref:`isochrone_objects`.  Each filter is
identified by a unique string, and an array of such strings
are passed into the Isochrone call. 

For example::
  
    # Use the HST WFC3-IR F127M and F153M filters, along with NIRC2 Kp
    filt_list = ['wfc3,ir,f127m', 'wfc3,ir,f153m', 'nirc2,Kp']
    my_iso = synthetic.IsochronePhot(logAge, AKs, dist, metallicity=0,
                            evo_model=evo_model, atm_func=atm_func,
                            red_law=red_law, filters=filt_list)
    
These strings follow the format ``<telescope/filter_set>,<filter>``.
Note that there is no space after the comma, and case matters.

Filter transmissions are stored in specific sub-directories
in ``SPISEA/filt_func``. The functions in ``spisea/filters.py`` tell
the :ref:`isochrone_objects` how to interpret the filter strings
and pull the corresponding transmission functions from the filt_func
directories.

Available filters:

* 2MASS
* CTIO_OSIRIS
* DeCam
* GAIA
* HAWK-I
* Hubble Space Telescope
* Johnson-Cousins
* Johnson-Glass
* JWST
* Keck NIRC
* Keck NIRC2
* NACO 
* PanStarrs 1
* UKIRT
* VISTA
* ZTF

  
Filter Sets
------------

   
**2MASS**

`Two-Micron Sky Survey <https://old.ipac.caltech.edu/2mass/>`_ Filters: J, H, Ks

Example: ``'2mass,H'``


**CTIO_OSIRIS**

`OSIRIS imager
<http://www.ctio.noao.edu/soar/content/ohio-state-infrared-imagerspectrograph-osiris>`_
on the CTIO telescope. Filters: H, K

Example: ``'ctio_osiris,H'``


**DeCam**

`Dark Energy Camera <http://www.ctio.noao.edu/noao/content/DECam-filter-information>`_
Filters: u, g, r, i, z, Y

Example: ``'decam,r'``

**GAIA**

The `GAIA Space Telescope filters <https://www.cosmos.esa.int/web/gaia/iow_20180316>`_.
Note that three sets are available: the pre-launch passbands used in DR1
(`Jordi+10
<https://ui.adsabs.harvard.edu/abs/2010A%26A...523A..48J/abstract>`_),
the passbands used for the DR2 published photometry, and
the *revised* DR2 passbands based on the DR2 data (October 2017).
ONLY THE REVISED DR2 PASSBANDS ARE SUPPORTED BY SPISEA.

Filters: G, Gbp, Grp

Example (gaia G filter from revised DR2 passbands):
``'gaia,dr2_rev,G'``

**HAWK-I**

The `High Acuity Wide Field K-band Imager
<https://www.eso.org/sci/facilities/paranal/instruments/hawki.html>`_
located on the ESO VLT (Unit Telescope 4).

Filters: J, H, Ks

Example: ``hawki,J``

**Hubble Space Telescope**

HST filters are defined by their `pysynphot OBSMODE strings <https://pysynphot.readthedocs.io/en/latest/appendixb.html#pysynphot-appendixb>`_. 

Example: ``'wfc3,ir,f125w'``


**Johnson-Cousins**

Johnson-Cousin filters (downloaded from
http://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves). Filters:
U, B, V, R, I

Example: ``'ubv,B'``

**Johnson-Glass**

Johnson-Glass passbands taken from `Bessell et al. 1988 <https://ui.adsabs.harvard.edu//#abs/1988PASP..100.1134B/abstract>`_
Filters: J, H, K

Example: ``'jg,K'``

**JWST**

JWST NIRCam filters, downloaded from `NIRCam website <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Filters#NIRCamFilters-filt_trans>`_. The filter functions in the nircam_throughputs/modAB_mean/nrc_plus_ote folder is used.

Filters: F070W, F090W,  F115W, F140M, F150W, F150W2, F162M, F164N, F182M, F187N, F200W, F210M, F212N, F250M, F277W, F300M, F322W2, F323N, F335M, F356W, F360M, F405N, F410M, F430M,  F444W, F460M, F466N, F470N, F480M, 

Example: ``'jwst,F356W'``
						

**Keck NIRC**

`NIRC1 filters <https://www2.keck.hawaii.edu/inst/nirc/>`_ on the Keck Telescope
Filters: H, K

Example: ``'nirc1,H'``


**Keck NIRC2**

`NIRC2 filters <https://www2.keck.hawaii.edu/inst/nirc2/filters.html>`_
on the Keck Telescope ()
Filters: J, H, Hcont, K, Kp, Ks, Kcont, Lp, Ms, Brgamma, FeII

Example: ``'nirc2,Ks'``


**NACO**

`ESO NACO filters <https://www.eso.org/sci/facilities/paranal/instruments/naco/inst/filters.html>`_
Filters: J, H, Ks

Example: ``'naco,H'``


**PanStarrs1**

PanStarrs 1 filters from `Tonry et al. 2012 <https://ui.adsabs.harvard.edu/#abs/arXiv:1203.0297>`_
Filters: g, r, i, z, y

Example: ``'ps1, g'``


**UKIRT**

UKIRT Telescope (http://www.ukidss.org/technical/photom/photom.html)
Filters: J, H, K

Example: ``'ukirt,K'``


**VISTA**

`VISTA Telescope <http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/filter-set>`_
Filters: Z, Y, J, H, K

Example: ``'vista,Y'``

**ZTF**

`ZTF Telescope <https://www.ztf.caltech.edu/page/technical>`_
Filters: g, r, i

Example: ``'ztf,g'``
