When calling synthetic.IsochronePhot, the user can specify the filters used for synthetic photometry. Each filter is identified by a unique string, and an array of such strings are passed into the <filters> variable (see Quick_Start_Make_Cluster.ipynb for example). These strings follow the format "<telescope/filter_set>,<filter>." Note that there is no space after the comma, and case matters.

Available filters:

2MASS
CTIO_OSIRIS
DeCam
Hubble Space Telescope
Johnson-Cousins
Johnson-Glass
JWST
Keck NIRC
Keck NIRC2
NACO 
PanStarrs 1
UKIRT
VISTA
GAIA

These filter functions are stored in the PopStar/filt_func directory. 

==========================
2MASS
==========================
Two-Micron Sky Survey (https://old.ipac.caltech.edu/2mass/)
Filters: J, H, Ks

Example: '2mass,H'

==========================
CTIO_OSIRIS
==========================
OSIRIS imager on the CTIO telescope (http://www.ctio.noao.edu/soar/content/ohio-state-infrared-imagerspectrograph-osiris)
Filters: H, K

Example: 'ctio_osiris,H'

==========================
DeCam
==========================
Dark Energy Camera (http://www.ctio.noao.edu/noao/content/DECam-filter-information)
Filters: u, g, r, i, z, Y

Example: 'decam,r'

==========================
Hubble Space Telescope
==========================
HST filters are defined by their pysynphot OBSMODE strings (https://pysynphot.readthedocs.io/en/latest/appendixb.html#pysynphot-appendixb). Note that the pysynphot cdbs/comp directory must be downloaded in order for this to work.

Example: WFC3-IR F125W filter: 'wfc3,ir,f125w'


==========================
Johnson-Cousins
==========================
Johnson-Cousin filters (downloaded from http://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves)
Filters: U, B, V, R, I

Example: 'ubv,B'

==========================
Johnson-Glass
==========================
Johnson-Glass passbands taken from Bessell et al. 1988 (https://ui.adsabs.harvard.edu//#abs/1988PASP..100.1134B/abstract)
Filters: J, H, K

Example: 'jg,K'

==========================
JWST
==========================
JWST NIRCam filters, downloaded from NIRCam website (https://jwst-docs.stsci.edu/display/JTI/NIRCam+Filters#NIRCamFilters-filt_trans). The filter functions in the nircam_throughputs/modAB_mean/nrc_plus_ote folder is used.

Filters: F070W, F090W,  F115W, F140M, F150W, F150W2, F162M, F164N, F182M, F187N, F200W, F210M, F212N, F250M, F277W, F300M, F322W2, F323N, F335M, F356W, F360M, F405N, F410M, F430M,  F444W, F460M, F466N, F470N, F480M, 

Example: 'jwst,F356W'
						
==========================
Keck NIRC
==========================
NIRC1 filters on the Keck Telescope (https://www2.keck.hawaii.edu/inst/nirc/)
Filters: H, K

Example: 'nirc1,H'

==========================
Keck NIRC2
==========================
NIRC2 filters on the Keck Telescope (https://www2.keck.hawaii.edu/inst/nirc2/filters.html)
Filters: J, H, Hcont, K, Kp, Ks, Kcont, Lp, Ms, Brgamma, FeII

Example: 'nirc2,Ks'

==========================
NACO
==========================
ESO NACO filters (https://www.eso.org/sci/facilities/paranal/instruments/naco/inst/filters.html)
Filters: J, H, K

Example: 'naco,H'

==========================
PanStarrs1
==========================
PanStarrs 1 filters from Tonry et al. 2012 (https://ui.adsabs.harvard.edu/#abs/arXiv:1203.0297)
Filters: g, r, i, z, y

Example: 'ps1, g'

==========================
UKIRT
==========================
UKIRT Telescope (http://www.ukidss.org/technical/photom/photom.html)
Filters: J, H, K

Example: 'ukirt,K'

==========================
VISTA
==========================
VISTA Telescope (http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/filter-set)
Filters: Z, Y, J, H, K

Example: "vista,Y"

==========================
GAIA
==========================
The GAIA Space Telescope filters
(https://www.cosmos.esa.int/web/gaia/iow_20180316).
Note that three sets are available: the pre-launch passbands used in DR1
(Jordi+10), the passbands used for the DR2 published photometry, and
the *revised* DR2 passbands based on the DR2 data (October 2017). The
user specifies which one they want by 'dr1', 'dr2', or 'dr2_rev', respectively.

To calculate synthetic fluxes, the dr2_rev passbands are advised.

Filters: G, Gbp, Grp

Example: 'gaia,dr2_rev,G'
