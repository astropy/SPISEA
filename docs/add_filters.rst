.. _add_filters:

========================================
Adding Photometric Filters
========================================
 If the user wants to add new photometric filters to PyPopStar, there are 3 main steps:

  1) Save the filter transmissions as text files in a new sub-directory in the top-level ``filt_func`` directory.
  2) Define a new function in ``filters.py`` that reads a unique filter string the user assigns to the new filters,
     and then read the appropriate data files from the new ``filt_func`` sub-directory. 
  3) Edit the ``get_filter_info()`` function  in ``synthetic.py`` to call the new function in ``filters.py`` when the new filter string is called.


Additional documentation on this is coming soon. In the meantime, let us know on the  Github `issue tracker
<https://github.com/astropy/PyPopStar/issues>`_ if you'd like to
implement new photometric filters.
