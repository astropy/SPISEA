---
title: 'PyPopStar: A Python-Based Simple Stellar Population Synthesis Code for Star Clusters'
tags:
  - Python
  - star clusters
  - stellar populations
authors:
  - name: Matthew W. Hosek Jr.
    orcid: 0000-0003-2874-1196
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Jessica R. Lu
    orcid: 0000-0001-9611-0009
    affiliation: 2
  - name: Casey Y. Lam
    orcid: 0000-0002-6406-1924
    affiliation: 2
  - name: Abhimat K. Gautam
    orcid: 0000-0002-2836-117X
    affiliation: 1
  - name: Kelly E. Lockhart
    orcid: 0000-0002-8130-1440
    affiliation: 3
  - name: Dongwon Kim
    orcid: 0000-0002-6658-5908
    affiliation: 2
  - name: Siyao Jia
    orcid: 0000-0001-5341-0765
    affiliation: 2
affiliations:
 - name: University of California, Los Angeles, Department of Astronomy, Los Angeles, CA 90095
   index: 1
 - name: University of California, Berkeley, Department of Astronomy, Berkeley, CA 94720
 index: 2
  - name: Harvard-Smithsonian Center for Astrophysics, 60 Garden Street, Cambridge, MA 02138, USA
   index: 3
date: 30 April 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: The Astronomical Journal <- The name of the AAS journal.
---

# Summary

The ability to simulate stellar populations is an essential tool for
interpreting observations of star clusters.
There are a range of codes
available for this task, each with a different set
of advantages.
Examples such as Starburst99 [@Leitherer:1999xy; @Vazquez:2005nr], 
PEGASE [@Fioc:1997ul; @Fioc:2019ix],
GALAXEV [@2003MNRAS.344.1000B], 
GALICS [@Hatton:2003rp], 
Flexible Stellar Population Synthesis (FSPS) [@Conroy:2009le; @Conroy:2010cr],
CIGALE [@Noll:2009if],
and Stochastically Lighting Up Galaxies (SLUG) [@da-Silva:2012hb; @Krumholz:2015dk]
specialize in modeling stellar populations on galactic scales,
with features such as models for photoionization and dust
attenuation and prescriptions for supernovae yields and 
chemical evolution.
Others such as
Binary Population and Spectral Synthesis (BPASS) [@Eldridge:2017rc; @Stanway:2018dn] and 
SYCLIST [@Georgy:2014wj]
offer advance treatment of binary stellar evolution and stellar
rotation.

However, a common drawback of these codes is that they limit the
user to fixed set of options when choosing the "ingredients" of the
stellar population, such as the Initial Mass Function (IMF),
extinction law, stellar multiplicity properties, and/or the
Initial-Final Mass Relation (IFMR).
This hinders the ability to model observations of star
clusters with enough flexibility to study these properties.

We present PyPopStar, an open-source Python
package that simulates simple stellar populations (single age and
metallicity).
The strength of PyPopStar is its modular interface which offers the
user control of 13 input properties, including (but not limited to) the
IMF, stellar multiplicity, extinction law, IMFR, and which
stellar evolution and atmosphere model grids are used.
These properties are constructed as code objects that are
straight-forward to manipulate, and allows the user to
create new sub-classes in order to implement new models and/or
functionalities that can be integrated with the rest of the code.

PyPopStar can be downloaded through
[GitHub](https://github.com/astropy/PyPopStar), and has documentation via [ReadtheDocs](https://pypopstar.readthedocs.io/en/latest/). 

# Features

PyPopStar generates single-age, single-metallicity populations (i.e. star clusters). It gives the user control over many parameters:

- Cluster characteristics (age, metallicity, mass, distance)
- Total extinction, differential extinction, and extinction law
- Stellar evolution and atmosphere models
- Stellar multiplicity and Initial Mass Function
- Initial-Final Mass Relation
- Photometric filters

Some of the tasks that PyPopStar can perform include:

- make a cluster isochrone in many filters using different stellar models
- make a star cluster at any age with an unusual IMF and unresolved multiplicity
- make a spectrum of a star cluster in integrated light

A more complete description of the code structure, underlying
calculations, and available features is provided in Hosek et
al. (2020, submitted to AJ).
In addition to the [documentation](https://pypopstar.readthedocs.io/en/latest/),
we also provide jupyter notebooks with a [quick-start tutorial](https://github.com/astropy/PyPopStar/blob/master/docs/Quick_Start_Make_Cluster.ipynb) and additional [user examples](https://github.com/astropy/PyPopStar/tree/master/docs/paper_examples). 

# Research

The range of use cases for PyPopStar is demonstrated in several
published studies: modeling the IMF of star clusters
[@Lu:2013wo; @Hosek:2019kk], measuring the extinction law in highly
reddened regions [@Hosek:2018lr], predicting black hole microlensing
rates [@Lam:2020vl], and calculating photometric transformations
between different filters at high extinction and with a non-standard
extinction law [@Krishna-Gautam:2019sj; @Chen:2019dp].
Additional code development is underway to expand PyPopStar's capabilities for additional science applications, and we encourage contributions and/or feature requests from the community.

# Acknowledgements


# References







