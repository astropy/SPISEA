# PyPopStar

PyPopStar is an python package that generates single-age, single-metallicity
populations (i.e. star clusters). It gives the user control over many parameters:

* Cluster characteristics (age, metallicity, mass, distance)
* Total extinction, differential extinction, and extinction law
* Stellar evolution and atmosphere models
* Stellar multiplicity and Initial Mass Function
* Initial-Final Mass Relation
* Photometric filters

Here is a brief list of things that PyPopStar can do:

* make a cluster isochrone in many filters using different stellar models
  
* make a star cluster at any age with an unusual IMF and unresolved multiplicity
  
* make a spectrum of a star cluster in integrated light

See [documentation](https://pypopstar.readthedocs.io/en/latest/) for details on 
[installing](https://pypopstar.readthedocs.io/en/latest/getting_started.html)
and running PyPopStar. We also provide jupyter notebooks with a 
[quick-start tutorial](https://github.com/astropy/PyPopStar/blob/master/docs/Quick_Start_Make_Cluster.ipynb)
and [additional examples](https://github.com/astropy/PyPopStar/tree/master/docs/paper_examples)
demonstrating how to use PyPopStar. 

## Contributions
We encourage contributions to PyPopStar, particular those that add support for star formation histories, new models, higher spectral resolution, etc. For feature additions, we ask that users fork or branch off of the development repository, make their changes, and then submit merge and pull requests.

## Other Resources

* _Astropy: http://www.astropy.org/
* _git: http://git-scm.com/
* _github: http://github.com
* _Cython: http://cython.org/

## License 
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
