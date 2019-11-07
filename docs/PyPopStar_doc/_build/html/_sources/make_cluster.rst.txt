.. _cluster_objects:

===============
Cluster Object
===============
The cluster classes are defined in popstar/synthetic.py. The primary
inputs to a cluster object is the cluster mass,
:ref:`isochrone_objects`, and :ref:`imf_objects`.  In addition, an
:ref:`ifmr_objects` may be defined to produce compact stellar remnants.

To make a cluster::
  
  from popstar import synthetic
  clust = synthetic.<Cluster_class>

See `Quick Start Example
<https://github.com/astropy/PyPopStar/blob/new_doc/docs/Quick_Start_Make_Cluster.ipynb>`_
for a detailed example for how to make a cluster with the different
sub-classes and interact with the subsequent output.

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



