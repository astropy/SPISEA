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
date: 13 August 2017
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: The Astronomical Journal <- The name of the AAS journal.
---

# Summary

-top level summary and statement of need here

-contributions, bugs/feature requests: github page

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

A more complete description of the code structure, underlying calculations, and available features is provided in Hosek et al. (2020, submitted to AJ).
Additional details can be found in the [documentation](https://pypopstar.readthedocs.io/en/latest/).
We also provide jupyter notebooks with a [quick-start tutorial](https://github.com/astropy/PyPopStar/blob/master/docs/Quick_Start_Make_Cluster.ipynb) and additional [user examples](https://github.com/astropy/PyPopStar/tree/master/docs/paper_examples). 

# Research

The range of use cases for PyPopStar is demonstrated in several published studies: modeling the IMF of star clusters \citep{Lu:2013wo, Hosek:2019kk}, measuring the extinction law in highly reddened regions \citep{Hosek:2018lr}, predicting black hole microlensing rates \citep{Lam:2020vl}, and calculating photometric transformations between different filters at high extinction and with a non-standard extinction law \citep{Krishna-Gautam:2019sj, Chen:2019dp}.
Additional code development is underway to expand PyPopStar's capabilities for additional science applications, and we encourage contributions and/or feature requests from the community.

# Acknowledgements


# References



# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Fenced code blocks are rendered with syntax highlighting:
```python
for n in range(10):
    yield f(n)
```	






