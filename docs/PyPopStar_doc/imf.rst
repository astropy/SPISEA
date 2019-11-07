.. _imf_objects:

==========
IMF Object
==========
Initial Mass Functions (IMFs) are defined as classes in
popstar/imf/imf.py. These can be defined by::

  from popstar.imf import imf
  imf_obj = imf.<class_name>

The IMF object is an input for the :ref:`cluster_objects`, and will
be used to draw the inital stellar mass distribution for the cluster.
A :ref:`multi_obj` may be passed to the IMF object to
form multiple star systems.

The IMFs are defined as broken power laws. Additional functional forms
may be implemented as sub-classes off of the base IMF class.


Base IMF Class
--------------
.. autoclass:: imf.imf.IMF
	       :members: generate_cluster, calc_multi

Broken Power-Law IMFs
---------------------
.. autoclass:: imf.imf.IMF_broken_powerlaw
	       :show-inheritance:
		  
.. autoclass:: imf.imf.IMFSalpeter1955
	       :show-inheritance:
		  
.. autoclass:: imf.imf.Miller_Scalo_1979
	       :show-inheritance:
   
.. autoclass:: imf.imf.Kennicutt_1983
	       :show-inheritance:
		  
.. autoclass:: imf.imf.Kroupa_2001
	       :show-inheritance:
		  
.. autoclass:: imf.imf.Weidner_Kroupa_2004
	       :show-inheritance:
