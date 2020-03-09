.. _cluster_objects:

===============
Cluster Object
===============
The cluster classes are defined in popstar/synthetic.py. The primary
inputs to a cluster object are the cluster mass,
:ref:`isochrone_objects`, and :ref:`imf_objects`.  In addition, an
:ref:`ifmr_objects` may be defined to produce compact stellar remnants.

An example of making a ResolvedCluster, assuming the isochrone object
has already been created::
 
  from popstar import synthetic, ifmr
  from popstar.imf import imf, multiplicity
  import numpy as np

  # Define stellar multiplicity properties. Here we 
  # use the default multiplicity object.
  # This is an input for the IMF object.
  imf_multi = multiplicity.MultiplicityUnresolved()

  # Define the IFMR. Here we use the default
  # IFMR object. 
  # If no IFMR is desired, set this variable 
  # to None
  my_ifmr = ifmr.IFMR()

  # Define the IMF. Here we'll use a broken
  # power-law with the parameters from 
  # Kroupa et al. (2001, MNRAS, 322, 231),
  # using the multiplicity we defined previously
  massLimits = np.array([0.08, 0.5, 1, 120]) # mass segments
  powers = np.array([-1.3, -2.3, -2.3]) # power-law exponents 
  my_imf = imf.IMF_broken_powerlaw(massLimits, powers, 
  imf_multi)

  # Define the cluster mass
  mass = 10**5 # Units: solar masses

  # Make the cluster object. We will assume that the isochrone object
  # is already defined as my_iso.
  cluster = synthetic.ResolvedCluster(my_iso, my_imf, mass, 
  ifmr=my_ifmr)


See `Quick Start Example
<https://github.com/astropy/PyPopStar/blob/new_doc/docs/Quick_Start_Make_Cluster.ipynb>`_
for a detailed example for how to make different
cluster sub-classes and interact with the subsequent output.

Base Cluster Class
----------------------------
.. autoclass:: synthetic.Cluster


Cluster Sub-classes
-----------------------
.. autoclass:: synthetic.ResolvedCluster
	       :show-inheritance:
	       
.. autoclass:: synthetic.ResolvedClusterDiffRedden
	       :show-inheritance:
	
.. autoclass:: synthetic.UnresolvedCluster
	       :show-inheritance:



