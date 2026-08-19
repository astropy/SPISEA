.. _multi_obj:

===========================
Stellar Multiplicity Object
===========================
The properties of multiple systems in the stellar population is
defined by the stellar multiplicity object. The multiplicity classes
are defined in ``spisea/imf/multiplicity.py``.

To call a multiplicity class and wire it into a cluster::

  from spisea.imf import imf, multiplicity
  from spisea import synthetic

  multi = multiplicity.<class_name>(...)
  imf_obj = imf.Kroupa_2001(multiplicity=multi)
  cluster = synthetic.ResolvedCluster(iso, imf_obj, Mcl)

The multiplicity object provides the following functions used by the IMF:

* ``multiplicity_fraction(mass)``
* ``companion_star_fraction(mass)``
* ``random_q(x, mass=None)`` — pass ``mass`` for mass-dependent q
  (brown-dwarf vs stellar). ``random_q(x)`` with no mass keeps the
  historical stellar-only power law.
* ``random_companion_count(x, CSF, MF, mass=None, rng=None)`` —
  companion-count policy, including the binaries-only BD cap when
  ``mass`` is given.
* attributes ``companion_max``, ``CSF_max``, ``q_min``

The user can choose either
an unresolved or a resolved multiplicity object. If a resolved
multiplicity object is selected, then orbital parameters are
assigned to each companion star (e.g semi-major axis, eccentricity,
inclination). These values are added as additional columns in the ``companions``
table off of the cluster object.

Note that the synthetic photometry
returned in the ``star_systems`` table off the cluster object is the
same for both unresolved and resolved multiplicity classes: it
represents the combined photometry of all stars within a given system.

For most selected evolution models, the companions are evolved as single stars.
To evolve binaries with mass exchange (does not support higher order multiples), you should use one of the ``MultiplicityResolved`` classes
and the ``COSMIC`` evolution model. 
See the example jupyter notebook `Cluster_w_COSMIC.ipynb <https://github.com/MovingUniverseLab/spisea/blob/main/docs/Cluster_w_COSMIC.ipynb>`_ for an example.
Note that currently COSMIC due to being external evolution is significantly slower than the other evolution options.


Unresolved Multiplicity Classes
------------------------------------------
.. autoclass:: imf.multiplicity.MultiplicityUnresolved
	       :members: companion_star_fraction,
			 multiplicity_fraction, random_q

.. autoclass:: imf.multiplicity.MultiplicityPiecewisePowerLaw
	       :show-inheritance:
	       :members: multiplicity_fraction, companion_star_fraction

.. autoclass:: imf.multiplicity.MultiplicityLogistic
	       :show-inheritance:
	       :members: multiplicity_fraction, companion_star_fraction

.. autoclass:: imf.multiplicity.MultiplicityUnresolvedOffner2023
	       :show-inheritance:
	       :members: multiplicity_fraction, companion_star_fraction,
			 q_power_at_mass, random_q, log_a_mean, a_mean,
			 sigma_log_a


Offner et al. 2023 multiplicity
------------------------------------------
The recommended multiplicity class to use is that derived from 
data summarized in Offner et al. (2023) (`arXiv:2203.10066
<https://arxiv.org/abs/2203.10066>`_; ADS
`2023ASPC..534..275O
<https://ui.adsabs.harvard.edu/abs/2023ASPC..534..275O>`_). Table 1
data: Zenodo `10.5281/zenodo.6628915
<https://doi.org/10.5281/zenodo.6628915>`_.
This class is not the default (for backwards compatability) but 
is strongly preferred.


The SPISEA v2.5 :class:`~imf.multiplicity.MultiplicityUnresolved` /
:class:`~imf.multiplicity.MultiplicityResolvedDK` objects remain the
default; but has known limitations in the brown dwarf regime.

Unresolved (companions, no orbits)::

  from spisea.imf import imf, multiplicity
  from spisea import synthetic

  multi = multiplicity.MultiplicityUnresolvedOffner2023()
  imf_obj = imf.Kroupa_2001(multiplicity=multi)
  cluster = synthetic.ResolvedCluster(iso, imf_obj, Mcl)

Resolved (same MF/CSF/q, plus mass-dependent separations)::

  multi = multiplicity.MultiplicityResolvedOffner2023()

The generic helpers :class:`~imf.multiplicity.MultiplicityLogistic` and
:class:`~imf.multiplicity.MultiplicityPiecewisePowerLaw` are available
for other surveys. Offner does **not** evaluate MF/CSF, :math:`\gamma`,
or separations as a piecewise interpolation of Table 1/2 knots.

Multiplicity and companion-star fractions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MF and CSF are a 4-parameter logistic in log-mass, fitted with equal
weight to the geom-mean MF/CF columns of Table 1:

.. math::

   f(M) = A + \frac{B - A}{1 + (M / M_0)^{-k}}

with :math:`(A, B, M_0, k) = (0.14, 0.99, 1.41, 1.25)` for MF and
:math:`(0.12, 2.35, 3.57, 0.96)` for CSF. The curve is C-infinity
smooth and saturates at :math:`B \approx 1` for MF. MF is clipped to
:math:`[0, 1]`. CSF is clipped to :math:`[0, \mathrm{CSF_{max}}]`,
raised to at least MF, and forced equal to MF for
:math:`M \le 0.08\,M_\odot` (binaries only).

SPISEA v2.5 uses :math:`\mathrm{MF} = 0.44\,M^{0.51}` (clipped
to 1) for arrays, plus a scalar-only brown-dwarf staircase
(0 / 8% / 16%). Cluster generation on that class still uses the
stellar power law for brown-dwarf primaries.

.. figure:: figures/mf_offner_vs_spisea2.5.png
   :alt: Offner 2023 vs SPISEA v2.5: multiplicity fraction vs primary mass
   :align: center

   Offner 2023 vs SPISEA v2.5: multiplicity fraction vs primary mass.
   Left: brown-dwarf zoom. Right: BD through early B. Solid: Offner
   logistic in log-mass. Dashed: SPISEA v2.5 :math:`0.44\,M^{0.51}`.
   Dotted: SPISEA v2.5 scalar BD staircase. Blue points: Offner
   Table 1. Shaded: :math:`M \le 0.08\,M_\odot`.

Mass-ratio index
~~~~~~~~~~~~~~~~
Companion mass ratios follow :math:`P(q) \propto q^{\gamma}` on
:math:`q_{\min} \le q \le 1` (default :math:`q_{\min} = 0.01`).
:math:`\gamma(M)` is an error-weighted logistic in log-mass fitted to
Table 1 :math:`\gamma_{\mathrm{trunc}}` (1–100 au):

.. math::

   \gamma(M) = A + \frac{B - A}{1 + (M / M_0)^{-k}}

with :math:`(A, B, M_0, k) = (6.6, -1.77, 0.0651, 0.629)`. Call
``q_power_at_mass(mass)`` or ``random_q(x, mass=...)``. Without
``mass``, ``random_q(x)`` keeps the historical stellar-only power law.

The err-weighted fit undershoots Fontanive et al. (2018)
:math:`8\pm 6\%` MF and :math:`\gamma = 4.8\pm 2.2`: at
:math:`0.033\,M_\odot`, :math:`\gamma \approx 3.3`. That is the
fit, not a bug. SPISEA v2.5 is a step: :math:`\gamma = 6.1` for
:math:`M \le 0.08\,M_\odot` (Fontanive) and :math:`-0.4` above.

.. figure:: figures/q_offner_vs_spisea2.5.png
   :alt: Offner 2023 vs SPISEA v2.5: mass-ratio index vs primary mass
   :align: center

   Offner 2023 vs SPISEA v2.5: mass-ratio index :math:`\gamma` vs
   primary mass. Solid: Offner error-weighted logistic. Dashed:
   SPISEA v2.5 step (6.1 below :math:`0.08\,M_\odot`, :math:`-0.4`
   above). Blue points: Table 1 :math:`\gamma_{\mathrm{trunc}}`.

.. figure:: figures/meanq_offner_vs_spisea2.5.png
   :alt: Offner 2023 vs SPISEA v2.5: mean mass ratio vs primary mass
   :align: center

   Offner 2023 vs SPISEA v2.5: mean mass ratio
   :math:`\langle q \rangle` on :math:`[0.01, 1]` implied by
   :math:`P(q)\propto q^{\gamma}`. Offner brown-dwarf companions
   are more equal-mass than SPISEA v2.5 stellar :math:`q`.

Characteristic separation
~~~~~~~~~~~~~~~~~~~~~~~~~
The characteristic :math:`\mu(a)` is a smooth broken power law in
:math:`\log_{10} a` vs :math:`\log_{10} M`, FGK-pulled, with smoothing
scale :math:`s = 0.1` dex. It is C-infinity (stable
:math:`\log\cosh`; not :math:`\log(\cosh x)` and not a hard
``where`` break):

.. math::

   v = \log_{10}(M / M_p),\quad
   y_p = \log_{10}(\mu_p)

.. math::

   \log_{10} a = y_p + \tfrac{1}{2}(\alpha_L+\alpha_R)\,v
     + \tfrac{1}{2}(\alpha_R-\alpha_L)\,s\,\log\cosh(v/s)

with :math:`\mu_p = 44.46` AU, :math:`M_p = 0.819\,M_\odot`,
:math:`\alpha_L = 1.005`, :math:`\alpha_R = -0.308`, :math:`s = 0.10`.
Linear-space :math:`a` is clipped above 0.1 AU. The implementation
uses the stable form
:math:`\log\cosh x = |x| + \log(1+e^{-2|x|}) - \log 2`.
``log_a_mean(mass)`` returns :math:`\log_{10}(a/\mathrm{AU})`;
``a_mean(mass)`` returns :math:`a` in AU.

The SPISEA v2.5 :class:`~imf.multiplicity.MultiplicityResolvedDK`
uses a Duchêne & Kraus (2013) broken power law in :math:`a` with a
brown-dwarf blend. That law is not meant for the BD regime.
``MultiplicityResolvedDK.a_mean`` / ``log_a_mean`` and
``sigma_log_a`` return the same characteristic mean and width that
``log_semimajoraxis`` draws from.

.. figure:: figures/sep_offner_vs_spisea2.5.png
   :alt: Offner 2023 vs SPISEA v2.5: characteristic separation vs primary mass
   :align: center

   Offner 2023 vs SPISEA v2.5: characteristic :math:`\mu(a)` vs
   primary mass. Solid: Offner smooth broken power law. Dashed:
   SPISEA v2.5 mean :math:`a`. Open circles: Table 1
   :math:`\tilde{a}_{\mathrm{all}}`. Filled squares: Table 2
   :math:`\mu` knots (4, 25, 40 AU). Offner BD binaries peak at a
   few AU.

Separation scatter
~~~~~~~~~~~~~~~~~~
:math:`\sigma(\log_{10} a)` is a 2-parameter logistic pinned at
0.7 / 1.5:

.. math::

   \sigma(M) = 0.7 + \frac{0.8}{1 + (M / 0.354)^{-6.05}}

i.e. :math:`(A, B, M_0, k) = (0.7, 1.5, 0.354, 6.05)`, clipped to
:math:`\ge 0.1`. Call ``sigma_log_a(mass)``. SPISEA v2.5 DK is a
linear fit in :math:`\log M` that saturates above
:math:`2.9\,M_\odot`; the dip near :math:`0.08\,M_\odot` is the
BD/stellar blend.

.. figure:: figures/sig_loga_offner_vs_spisea2.5.png
   :alt: Offner 2023 vs SPISEA v2.5: sigma of log10 a vs primary mass
   :align: center

   Offner 2023 vs SPISEA v2.5: :math:`\sigma(\log_{10} a)` vs
   primary mass. Solid: Offner 2-parameter logistic. Dashed:
   SPISEA v2.5 DK. Filled squares: Table 2 knots (0.7, 1.3, 1.5).

Resolved draws
~~~~~~~~~~~~~~
:class:`~imf.multiplicity.MultiplicityResolvedOffner2023` draws
:math:`\log_{10}(a/\mathrm{AU})` from a truncated lognormal with
``loc = log_a_mean(mass)`` and ``scale = sigma_log_a(mass)``,
truncated to 0.01–2000 AU (same limits as
:class:`~imf.multiplicity.MultiplicityResolvedDK`). Eccentricity
and Keplerian angles still follow Duchêne & Kraus (2013)
(:math:`f(e)=2e`, random inclination and angles).

Brown-dwarf policy and Table 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Primaries at or below :math:`0.08\,M_\odot` are binaries only:
CSF = MF and the companion count is capped at 1.

Table 1 FGKM MF/CF exclude brown-dwarf companions
(:math:`M_{\mathrm{comp}} > 0.075\,M_\odot` for FGKM;
OBA use :math:`q > 0.1`). SPISEA still draws companions down to
``q_min`` (default 0.01), so some stellar primaries get BD
secondaries (~4% for solar-type). Do not read the simulated
stellar-primary MF as a stellar-companion-only statistic.

Reproducing the figures
~~~~~~~~~~~~~~~~~~~~~~~
The comparison figures are generated from the multiplicity object
methods (``multiplicity_fraction``, ``q_power_at_mass``,
``a_mean``, ``sigma_log_a``) so they cannot drift from the code.
From the repository root::

   python docs/figures/plot_mf_offner_vs_spisea2.5.py
   python docs/figures/plot_q_offner_vs_spisea2.5.py
   python docs/figures/plot_meanq_offner_vs_spisea2.5.py
   python docs/figures/plot_sep_offner_vs_spisea2.5.py
   python docs/figures/plot_sig_loga_offner_vs_spisea2.5.py

``python docs/figures/plot_q_sep_offner_vs_spisea2.5.py`` still writes
the last four PNGs by running those scripts.


Resolved Multiplicity Classes
------------------------------------------
.. autoclass:: imf.multiplicity.MultiplicityResolvedDK
	       :show-inheritance:
	       :members: log_semimajoraxis, log_a_mean, a_mean, sigma_log_a

.. autoclass:: imf.multiplicity.MultiplicityResolvedOffner2023
	       :show-inheritance:
	       :members: log_semimajoraxis, log_a_mean, a_mean, sigma_log_a

