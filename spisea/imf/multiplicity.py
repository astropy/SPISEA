import numpy as np
import astropy.modeling
from random import choice
from scipy.stats import truncnorm

defaultMF_amp = 0.44
defaultMF_power = 0.51
defaultCSF_amp = 0.50
defaultCSF_power = 0.45
defaultCSF_max = 3
defaultq_power = -0.4
defaultq_min = 0.01
default_aMean = 100.0 # log (AU)
default_aSigma = 0.1  # log (AU)

# Hydrogen-burning limit used for BD-primary (binaries-only) logic.
# Offner et al. 2023 use M_comp > 0.075 Msun as the MS companion cut;
# SPISEA keeps 0.08 Msun for consistency with existing BD handling.
H_BURNING_MASS = 0.08

# Fontanive et al. (2018) mass-ratio power-law index used for BD primaries
# in the Lu et al. (2013) MultiplicityUnresolved implementation.
FONTANIVE2018_BD_Q_POWER = 6.1

# Eventually we should add in separation properties. (a_mean, a_sigma)

# Equal-weight logistic-in-log-mass fit to Offner et al. 2023 Table 1
# geom-mean (M, MF) and (M, CF) points:
#   y(M) = A + (B - A) / (1 + (M / M0)**(-k))
OFFNER2023_MF_A = 0.14
OFFNER2023_MF_B = 0.99
OFFNER2023_MF_M0 = 1.41
OFFNER2023_MF_K = 1.25
OFFNER2023_CSF_A = 0.12
OFFNER2023_CSF_B = 2.35
OFFNER2023_CSF_M0 = 3.57
OFFNER2023_CSF_K = 0.96


# Error-weighted logistic in log-mass for Table 1 γ_trunc (1–100 au).
# The model is this logistic, not interpolation of the arrays below.
OFFNER2023_Q_A = 6.6
OFFNER2023_Q_B = -1.77
OFFNER2023_Q_M0 = 0.0651
OFFNER2023_Q_K = 0.629

# Smooth broken power law in log10(a) vs log10(M), s=0.1 dex, FGK-pulled.
OFFNER2023_A_MUP = 44.46
OFFNER2023_A_MP = 0.819
OFFNER2023_A_ALPHAL = 1.005
OFFNER2023_A_ALPHAR = -0.308
OFFNER2023_A_S = 0.10
OFFNER2023_A_MIN = 0.1

# 2-parameter logistic for σ(log10 a); floors/ceilings pinned at 0.7 / 1.5.
OFFNER2023_SIG_A = 0.7
OFFNER2023_SIG_B = 1.5
OFFNER2023_SIG_M0 = 0.354
OFFNER2023_SIG_K = 6.05


class _ResolvedOrbitalMixin(object):
    """Eccentricity and Keplerian angles shared by resolved multiplicity classes."""

    def random_e(self, x):
        """
        Generate random eccentricity from the inverse of the CDF where the PDF is f(e) = 2e from Duchene and Kraus 2013
        
        Parameters
        ----------
        x : float or array_like
            Random number between 0 and 1.

        Returns
        -------
        e : float or array_like
            eccentricity
        """
        e = np.sqrt(x)
        
        return e
    
    def random_keplarian_parameters(self, x, y, z):
        """
        Generate random incliniation and angles of binary system
        
        Parameters
        ----------
        x : float or array_like
            Random number between 0 and 1.
            
        y : float or array_like
            Random number between 0 and 1.
            
        z : float or array_like
            Random number between 0 and 1.

        Returns
        -------
        inclination : float or array_like
            Angle of inclination
                    
        Omega : float or array_like
            Big Omega: one other angle of the system
        
        omega : float or array_like
            Final angle of the system
        """
        sign = np.array([choice([-1,1]) for i in range(len(x))])
        x = sign*x
        inclination = np.arccos(x)*180/np.pi #inclination angle in degrees
        
        Omega = 360*y
        omega = 360*z
        
        return inclination, Omega, omega

class MultiplicityUnresolved(object):
    """
    The properties of stellar companions (see notes below). 
    The default parameters are as described in 
    `Lu et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013ApJ...764..155L/abstract>`_.
    These parameters are most appropriate for stellar populations
    with ages <10 Myr.

    Notes
    -----
    The number of stellar companions, their masses, and separations
    are be described by the following functions:

    **Multiplicity Fraction** -- the number of stellar systems that host 
    multiple stars. In other words, the number of primary stars with
    companions. The multiplicity fraction (MF) is typically described
    as::
                            B + T + Q + ...
                MF =     ---------------------
                          S + B + T + Q + ...

    where S = single, B = binary, T = triple, Q = quadruple, etc.
    The MF also changes with mass and this dependency can be 
    described as a power-law::
            
                MF(mass) = MF_amp * (mass ** MF_power)

    However, in the brown dwarf mass regime, it is currently recognized
    that only binaries are possible, and the MF decreases dissimilarly
    to higher masses (> 0.08 solar masses). The values for this range
    are given by Aberasturi et al. (2014) and Fontanive et al. (2023).

    **Companion Star Fraction** -- the expected number of companions in
    a multiple system. The companion star fraction (CSF) also 
    changes with mass and this dependency can be described as
    a power-law::
                
                CSF(mass) = CSF_amp * (mass ** CSF_power)

    The companion star fraction is clipped to some maximum
    value, CSF_max. The actual number of companions is drawn 
    from a Poisson distribution with an expectation value of CSF.

    In the brown dwarf regime we impose an assumption that only
    binary systems are possible due to current literature trends.

    **Mass Ratio (Q)** -- The ratio between the companion star
    mass and primary star mass, Q = (m_comp / m_prim ) has
    a probability density function described by a powerlaw::

                P(Q) = Q ** q_power  for q_min <= Q <= 1

    Current observations show no significant mass dependence.
        
    Parameters
    ----------
    MF_amp : float, optional
        The amplitude of the power-law describing the Multiplicity 
        Fraction as a function of stellar mass. 

    MF_power : float, optional
        The power of the power-law describing the Multiplicity
        Fraction as a function of stellar mass.

    CSF_amp : float, optional
        The amplitude of the power-law describing the companion star 
        fraction as a function of stellar mass.

    CSF_power : float, optional
        The power of the power-law describing the companion star 
        fraction as a function of stellar mass.

    CSF_max : float, optional
        The maximum allowed companion star fraction, which is the
        expectation value for the number of companion stars. Given
        a CSF_max = 3, some systems will still have more than 3 
        companions.

    q_power : float, optional
        The power of the power-law describing the probability
        density function for the mass ratio.

    q_min : float, optional
        The minimum allowed Q value for the probability
        density function of the mass ratio.
    
    companion_max : bool, optional
        Sets CSF_max is the max as the max number of companions.
        Default False.

    binary_only_mass_max : float, optional
        Primary mass in solar masses (Msun) at and below which systems
        are restricted to at most one companion (CSF = MF). Default is
        0.08 Msun.
    
    """
    def __init__(self, 
                 MF_amp=0.44, MF_power=0.51,
                 CSF_amp=0.50, CSF_power=0.45, CSF_max=3,
                 q_power=-0.4, q_min=0.01, companion_max = False,
                 binary_only_mass_max=H_BURNING_MASS):
         
        self.MF_amp = MF_amp
        self.MF_pow = MF_power
        self.CSF_amp = CSF_amp
        self.CSF_pow = CSF_power
        self.CSF_max = CSF_max
        self.q_pow = q_power
        self.q_min = q_min
        self.companion_max = companion_max
        self.binary_only_mass_max = binary_only_mass_max

    def multiplicity_fraction(self, mass):
        """
        Given a star's mass, determine the probability that the star is in a
        multiple system (multiplicity fraction = MF).

        Modified to allow binary fraction to decrease in brown dwarf regime.
        Supported by Aberasturi et al. (2014) and Fontanive et al. (2018).

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        mf : float or ndarray
            Multiplicity fraction, dimensionless, in [0, 1].
            The fraction of stars at this mass that will have one or
            more companions. Python float if ``mass`` is scalar,
            ndarray otherwise.
        """
        # Multiplicity Fraction
        mf = self.MF_amp * mass ** self.MF_pow

        if np.isscalar(mf):
            if mf > 1:
                mf = 1
            # physically override mf for brown dwarfs
            if (mass <= 0.08) & (mass > 0.06):
                mf = 0.16
            if (mass <= 0.06) & (mass > 0.02):
                mf = 0.08
            if (mass < 0.02):
                mf = 0
        else:
            mf[mf > 1] = 1

        return mf

    def companion_star_fraction(self, mass):
        """
        Given a star's mass, determine the average number of
        companion stars (companion star fraction = CSF). For
        brown dwarfs we impose a hard limit of one companion.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        csf : float or ndarray
            Companion star fraction, the expected number of companions
            for a star at this mass. Dimensionless mean companion count
            (not bounded by 1). Python float if ``mass`` is scalar,
            ndarray otherwise.
        """
        # Companion Star Fraction
        csf = self.CSF_amp * mass ** self.CSF_pow
        
        if np.isscalar(csf):
            if csf > self.CSF_max:
                csf = self.CSF_max
            if (mass <= self.binary_only_mass_max):
                csf = self.multiplicity_fraction(mass)
        else:
            csf[csf > self.CSF_max] = self.CSF_max
            bd = mass <= self.binary_only_mass_max
            csf[bd] = self.multiplicity_fraction(mass[bd])

        return csf

    def q_power_at_mass(self, mass):
        """
        Mass-ratio power-law index, P(q) ∝ q ** q_power.

        Lu et al. (2013) use a single ``q_power`` for stellar primaries.
        Brown-dwarf primaries (M <= binary_only_mass_max) use
        gamma = 6.1 from Fontanive et al. (2018), matching the
        companion-mass draw previously special-cased in ``imf.calc_multi``.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        q_power : float or ndarray
            Mass-ratio power-law index γ, dimensionless.
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
        q_pow = np.full(mass_arr.shape, self.q_pow, dtype=float)
        q_pow[mass_arr <= self.binary_only_mass_max] = FONTANIVE2018_BD_Q_POWER
        if np.isscalar(mass):
            return float(q_pow[0])
        return q_pow

    def random_q(self, x, mass=None):
        """
        Generative function for companion mass ratio, equivalent
        to the inverse of the CDF.

            `q = m_companion / m_primary`
            `P(q) = q ** beta`    for q_min <= q <= 1

        Parameters
        ----------
        x : float or array_like
            Uniform random draw, dimensionless, in [0, 1]. Inverse CDF
            sample for q.

        mass : float or array_like, optional
            Primary mass in solar masses (Msun). If given, the
            power-law index is ``q_power_at_mass(mass)`` (brown-dwarf
            vs stellar for the SPISEA v2.5 default; mass-dependent for
            Offner et al. 2023). If omitted, ``self.q_pow`` is used
            for all companions.

        Returns
        -------
        q : float or ndarray
            Companion mass ratio m_comp/m_prim, dimensionless, in
            [q_min, 1]. Python float if ``x`` is scalar, ndarray
            otherwise.
        """
        if mass is None:
            return _q_from_powerlaw(x, self.q_pow, self.q_min)
        return _q_from_powerlaw(x, self.q_power_at_mass(mass), self.q_min)

    def random_is_multiple(self, x, MF):
        """
        Helper function: determine if star is in multiple system.
        """
        return x < MF

    def random_companion_count(self, x, CSF, MF, mass=None, rng=None):
        """
        Number of companions for primaries already identified as multiple.

        The count is drawn from a Poisson with expectation CSF/MF - 1,
        then 1 is added so every multiple has at least one companion.
        ``x`` is unused and kept for API compatibility.

        Parameters
        ----------
        x : float or array_like
            Unused (historical signature). Dimensionless uniform
            draw in [0, 1] if provided.
        CSF : float or array_like
            Companion star fraction, dimensionless mean companion
            count (not bounded by 1).
        MF : float or array_like
            Multiplicity fraction, dimensionless, in [0, 1].
        mass : float or array_like, optional
            Primary mass in solar masses (Msun). If given, primaries
            at or below ``binary_only_mass_max`` are limited to one
            companion. Cluster generation always passes mass so
            subclasses can override the BD companion-count policy here.
        rng : numpy.random.Generator, optional
            Random generator. If omitted, uses ``numpy.random`` (the
            historical scalar helper).

        Returns
        -------
        n_comp : int or ndarray of int
            Number of companions, integer count. Python int if ``CSF``
            and ``MF`` are scalar, ndarray otherwise.
        """
        return_scalar = np.isscalar(CSF) and np.isscalar(MF)
        if return_scalar and rng is None:
            if MF <= 0:
                return 0
            n_comp = 1 + np.random.poisson((CSF / MF) - 1)
            if self.companion_max and n_comp > self.CSF_max:
                n_comp = self.CSF_max
            if mass is not None and np.asarray(mass, dtype=float).reshape(-1)[0] <= self.binary_only_mass_max:
                n_comp = min(int(n_comp), 1)
            return int(n_comp)

        CSF = np.atleast_1d(np.asarray(CSF, dtype=float))
        MF = np.atleast_1d(np.asarray(MF, dtype=float))
        if rng is None:
            n_comp = 1 + np.random.poisson((CSF / MF) - 1)
        else:
            n_comp = 1 + rng.poisson((CSF / MF) - 1)
        if self.companion_max:
            n_comp = np.minimum(n_comp, self.CSF_max)
        if mass is not None:
            mass = np.atleast_1d(np.asarray(mass, dtype=float))
            bd = mass <= self.binary_only_mass_max
            n_comp[bd & (n_comp > 1)] = 1
        if return_scalar:
            return int(n_comp[0])
        return n_comp

    def draw_n_companions(self, mass, CSF, MF, rng):
        """
        Vectorized companion counts for primaries that are already
        identified as multiple. Delegates to :meth:`random_companion_count`
        with ``mass`` so BD companion-count policy lives on this object.

        Parameters
        ----------
        mass : array_like
            Primary masses in solar masses (Msun) of systems already
            identified as multiple.
        CSF : array_like
            Companion star fraction at each primary, dimensionless
            mean companion count (not bounded by 1).
        MF : array_like
            Multiplicity fraction at each primary, dimensionless,
            in [0, 1].
        rng : numpy.random.Generator
            Random generator used for the Poisson companion-count draw.

        Returns
        -------
        n_comp : ndarray of int
            Number of companions per primary, integer count, shape
            matching ``mass``.
        """
        return np.atleast_1d(
            self.random_companion_count(None, CSF, MF, mass=mass, rng=rng))

    def _q_values_for_primaries(self, prim_subset, n_comp, rng):
        """
        Draw mass ratios for ``n_comp`` companions of each primary.

        The stellar / brown-dwarf split (two separate RNG draws) preserves
        the historical Lu et al. (2013) random sequence used by
        ``imf.calc_multi``.

        Parameters
        ----------
        prim_subset : array_like
            Primary masses in solar masses (Msun) for this companion-count
            group.
        n_comp : int
            Number of companions per primary, integer count.
        rng : numpy.random.Generator
            Random generator used for inverse-CDF q draws.

        Returns
        -------
        q_values : ndarray
            Companion mass ratios m_comp/m_prim, dimensionless, in
            [q_min, 1]. Shape (len(prim_subset), n_comp).
        """
        q_values = np.empty((len(prim_subset), n_comp))
        bd_mask = prim_subset <= self.binary_only_mass_max
        star_mask = ~bd_mask
        if np.any(star_mask):
            q_values[star_mask] = self.random_q(
                rng.random((star_mask.sum(), n_comp)))
        if np.any(bd_mask):
            q_values[bd_mask] = self.random_q(
                rng.random((bd_mask.sum(), n_comp)),
                mass=prim_subset[bd_mask])
        return q_values

    def draw_companion_masses(self, primary_masses, is_multiple, CSF, MF,
                              rng, mass_min):
        """
        Assign companion masses for a set of primaries.

        This is the multiplicity-object entry point used by
        ``IMF.calc_multi``. Companion-mass draws, including brown-dwarf
        q distributions and the binaries-only BD companion count, live
        here rather than in ``imf.py``.

        Parameters
        ----------
        primary_masses : array_like
            Primary masses in solar masses (Msun).
        is_multiple : array_like of bool
            True for primaries drawn as multiple systems.
        CSF : array_like
            Companion star fraction at each primary, dimensionless
            mean companion count (not bounded by 1).
        MF : array_like
            Multiplicity fraction at each primary, dimensionless,
            in [0, 1].
        rng : numpy.random.Generator
            Random generator.
        mass_min : float
            Minimum companion mass in solar masses (Msun); lighter
            companions are masked.

        Returns
        -------
        comp_masses : numpy.ma.MaskedArray
            Companion masses in solar masses (Msun), shape
            (n_primaries, max_n_comp).
        system_masses : ndarray
            Primary plus unmasked companion mass, in solar masses
            (Msun).
        is_multiple : ndarray of bool
            Updated multiplicity flags after masking sub-minimum
            companions.
        """
        primary_masses = np.asarray(primary_masses, dtype=float)
        is_multiple = np.asarray(is_multiple, dtype=bool)
        CSF = np.asarray(CSF, dtype=float)
        MF = np.asarray(MF, dtype=float)
        system_masses = primary_masses.copy()
        multiple_idx = np.where(is_multiple)[0]
        primary = primary_masses[multiple_idx]
        n_comp = self.draw_n_companions(
            primary, CSF[multiple_idx], MF[multiple_idx], rng)

        if len(multiple_idx) == 0:
            comp_masses = np.zeros((len(primary_masses), 1))
            comp_masses = np.ma.MaskedArray(comp_masses, mask=comp_masses < mass_min)
            return comp_masses, system_masses, is_multiple

        n_unique = np.unique(n_comp)
        n_indices = [np.where(n_comp == i)[0] for i in n_unique]
        comp_masses = np.zeros((len(primary_masses), int(np.max(n_unique))))

        for n_c, idx in zip(n_unique, n_indices):
            prim_subset = primary[idx]
            q_values = self._q_values_for_primaries(prim_subset, int(n_c), rng)
            m_comp = q_values * prim_subset[:, np.newaxis]
            comp_masses[multiple_idx[idx], :int(n_c)] = m_comp

        comp_masses = np.ma.MaskedArray(comp_masses, mask=comp_masses < mass_min)
        system_masses[multiple_idx] += comp_masses[multiple_idx].sum(axis=1)
        is_multiple = np.any(~comp_masses.mask, axis=1)
        return comp_masses, system_masses, is_multiple
    
class MultiplicityResolvedDK(MultiplicityUnresolved, _ResolvedOrbitalMixin):
    """
    Sub-class of MultiplicityUnresolved that adds semimajor axis and eccentricity information 
    for multiple objects from distributions described in Duchene and Kraus 2013

    For brown dwarf regime, mean separation and std are given by Fontanive et al. (2018).

    Parameters
    --------------
    a_amp: float, optional
        Ampltiude of the broken power law describing the log_semimajoraxis
        
    a_break: float, optional
        Break location on the x-axis of the broken power law describing the log_semimajoraxis
    
    a_slope1: float, optional
        Slope of the left side of the broken power law describing the log_semimajoraxis
        
    a_slope2: float, optional
        Slope of the right side of the broken power law describing the log_semimajoraxis
        
    a_std_slope: float, optional
        Slope of the line that fit sigma_log_semimajoraxis vs log_mass
    
    a_std_intercept: float, optional
        Intercept of the line that fit sigma_log_semimajoraxis vs log_mass        
    """
    def __init__(self, a_amp = 379.79953034, a_break = 4.90441533, a_slope1 = -1.80171539, 
                 a_slope2 = 4.23325571, a_std_slope = 1.19713084, a_std_intercept = 1.28974264, **kwargs):
        super(MultiplicityResolvedDK, self).__init__(**kwargs)
        self.a_amp = a_amp
        self.a_break = a_break
        self.a_slope1 = a_slope1
        self.a_slope2 = a_slope2
        self.a_std_slope = a_std_slope
        self.a_std_intercept = a_std_intercept
    
    def log_semimajoraxis(self, mass):
        """
        Generate the semimajor axis for a given mass. The mean and standard deviation of a given mass are determined 
        by fitting the data from fitting the semimajor axis data as a function of mass in table 1 of Duchene and Kraus 2013.
        Then a random semimajor axis is drawn from a log normal distribution with that mean and standard deviation.

        The brown dwarf range is covered by mass-dependent scaling of both the characteristic separation and dispersion
        matching trends described in Fontanive et al. (2018).

        Parameters
        ----------
        mass : array-like
            Mass array of primary star

        Returns
        -------
        log_semimajoraxis : array-like
            Log of the semimajor axis/separation between the stars in units of AU
        """
        mass = np.atleast_1d(mass)
        logm = np.log10(mass)

        # Stellar mean and std (Duchene & Kraus 2013)
        a_mean_func = astropy.modeling.powerlaws.BrokenPowerLaw1D(amplitude=self.a_amp, x_break=self.a_break,
                                                                  alpha_1=self.a_slope1, alpha_2=self.a_slope2)
        log_a_mean_star = np.log10(a_mean_func(mass))  # mean log(a)
        log_a_std_func = astropy.modeling.models.Linear1D(slope=self.a_std_slope, intercept=self.a_std_intercept)
        log_a_std_star = log_a_std_func(logm)  # sigma_log(a)
        log_a_std_star[mass >= 2.9] = log_a_std_func(np.log10(2.9))  # sigma_log(a)
        log_a_std_star = np.clip(log_a_std_star, 0.1, None)

        # BD mean and std (Fontanive+18): interpolated over substellar range
        log_a_mean_bd = np.interp(
            logm,
            [np.log10(0.01), np.log10(0.08)],
            [np.log10(2.5), np.log10(8.0)]
        )
        log_a_std_bd = np.interp(
            logm,
            [np.log10(0.01), np.log10(0.08)],
            [0.25, 0.5]
        )

        # Sigmoid blend: smoothly transitions from BD to stellar regime at 0.08 M_sun
        w = 1.0 / (1.0 + np.exp(-(logm - np.log10(0.08)) / 0.15))
        log_a_mean = (1 - w) * log_a_mean_bd + w * log_a_mean_star
        log_a_std = (1 - w) * log_a_std_bd + w * log_a_std_star

        # Trunc normal distribution between log10(0.01) AU and log10(2000) AU
        log_a_lower = np.log10(0.01)
        log_a_upper = np.log10(2000)
        a_lower_std = (log_a_lower - log_a_mean) / log_a_std
        a_upper_std = (log_a_upper - log_a_mean) / log_a_std

        log_semimajoraxis = truncnorm.rvs(a_lower_std, a_upper_std, loc=log_a_mean, scale=log_a_std)
        return log_semimajoraxis

class MultiplicityPiecewisePowerLaw(MultiplicityUnresolved):
    """
    Multiplicity described by a piecewise power law in primary mass.

    On each mass segment i, with edges ``mass_limits[i] <= M < mass_limits[i+1]``::

        MF(M)  = MF_amp[i]  * M ** MF_power[i]
        CSF(M) = CSF_amp[i] * M ** CSF_power[i]

    MF is clipped to [0, 1]. CSF is clipped to [0, CSF_max] and forced
    equal to MF for primaries at or below ``binary_only_mass_max``
    (binaries only). CSF is also raised to at least MF so the Poisson
    companion-count draw is well defined.

    Parameters
    ----------
    mass_limits : array_like
        Segment edges in solar masses (Msun), length N+1, strictly
        increasing.
    MF_amps : array_like
        Length-N amplitudes for the multiplicity fraction,
        dimensionless (units of MF / Msun**MF_power).
    MF_powers : array_like
        Length-N power-law indices for the multiplicity fraction,
        dimensionless.
    CSF_amps : array_like
        Length-N amplitudes for the companion star fraction,
        dimensionless (mean companion count / Msun**CSF_power).
    CSF_powers : array_like
        Length-N power-law indices for the companion star fraction,
        dimensionless.
    CSF_max : float, optional
        Maximum companion star fraction, dimensionless mean companion
        count (not bounded by 1). Default 3.
    q_power : float, optional
        Mass-ratio power-law index, dimensionless. Default -0.4.
    q_min : float, optional
        Minimum mass ratio m_comp/m_prim, dimensionless. Default 0.01.
    companion_max : bool, optional
        If True, cap companion counts at CSF_max. Default False.
    binary_only_mass_max : float, optional
        Primary mass in solar masses (Msun) at and below which systems
        are binaries only. Default 0.08 Msun.
    """
    def __init__(self, mass_limits, MF_amps, MF_powers, CSF_amps, CSF_powers,
                 CSF_max=3, q_power=-0.4, q_min=0.01, companion_max=False,
                 binary_only_mass_max=H_BURNING_MASS):
        mass_limits = np.asarray(mass_limits, dtype=float)
        MF_amps = np.asarray(MF_amps, dtype=float)
        MF_powers = np.asarray(MF_powers, dtype=float)
        CSF_amps = np.asarray(CSF_amps, dtype=float)
        CSF_powers = np.asarray(CSF_powers, dtype=float)
        nseg = len(MF_amps)
        if len(mass_limits) != nseg + 1:
            raise ValueError('len(mass_limits) must be len(MF_amps) + 1')
        if not (len(MF_powers) == len(CSF_amps) == len(CSF_powers) == nseg):
            raise ValueError('MF/CSF amplitude and power arrays must have equal length')
        if np.any(np.diff(mass_limits) <= 0):
            raise ValueError('mass_limits must be strictly increasing')

        super(MultiplicityPiecewisePowerLaw, self).__init__(
            MF_amp=MF_amps[-1], MF_power=MF_powers[-1],
            CSF_amp=CSF_amps[-1], CSF_power=CSF_powers[-1],
            CSF_max=CSF_max, q_power=q_power, q_min=q_min,
            companion_max=companion_max,
            binary_only_mass_max=binary_only_mass_max)
        self.mass_limits = mass_limits
        self.MF_amps = MF_amps
        self.MF_powers = MF_powers
        self.CSF_amps = CSF_amps
        self.CSF_powers = CSF_powers

    def multiplicity_fraction(self, mass):
        """
        Multiplicity fraction as a piecewise power law in primary mass.
        Clipped to [0, 1].

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        mf : float or ndarray
            Multiplicity fraction, dimensionless, in [0, 1].
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        return _piecewise_powerlaw(
            mass, self.mass_limits, self.MF_amps, self.MF_powers,
            clip_min=0.0, clip_max=1.0)

    def companion_star_fraction(self, mass):
        """
        Companion star fraction as a piecewise power law in primary mass.

        Clipped to [0, CSF_max], raised to at least MF, and set equal
        to MF for primaries at or below ``binary_only_mass_max``.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        csf : float or ndarray
            Companion star fraction, dimensionless mean companion
            count (not bounded by 1). Python float if ``mass`` is
            scalar, ndarray otherwise.
        """
        return_scalar = np.isscalar(mass)
        mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
        csf = _piecewise_powerlaw(
            mass_arr, self.mass_limits, self.CSF_amps, self.CSF_powers,
            clip_min=0.0, clip_max=self.CSF_max)
        mf = _piecewise_powerlaw(
            mass_arr, self.mass_limits, self.MF_amps, self.MF_powers,
            clip_min=0.0, clip_max=1.0)
        csf = np.maximum(csf, mf)
        bd = mass_arr <= self.binary_only_mass_max
        csf[bd] = mf[bd]
        if return_scalar:
            return float(csf[0])
        return csf


class MultiplicityLogistic(MultiplicityUnresolved):
    """
    Multiplicity described by a logistic in log primary mass.

        f(M) = A + (B - A) / (1 + (M / M0)**(-k))

    As M → 0, f → A; as M → ∞, f → B. The same functional form is
    used for MF and CSF with independent coefficients.
    ``MultiplicityUnresolvedOffner2023`` uses this class with
    coefficients fitted to Offner et al. (2023) Table 1.

    MF is clipped to [0, 1]. CSF is clipped to [0, CSF_max], raised
    to at least MF, and forced equal to MF for primaries at or below
    ``binary_only_mass_max`` (binaries only).

    Parameters
    ----------
    MF_A, MF_B : float
        Low-mass and high-mass asymptotes of the multiplicity-fraction
        logistic, dimensionless (MF).
    MF_M0 : float
        Characteristic mass of the MF logistic, in solar masses (Msun).
    MF_k : float
        MF logistic slope, dimensionless.
    CSF_A, CSF_B : float
        Low-mass and high-mass asymptotes of the companion-star-fraction
        logistic, dimensionless mean companion count (not bounded by 1).
    CSF_M0 : float
        Characteristic mass of the CSF logistic, in solar masses (Msun).
    CSF_k : float
        CSF logistic slope, dimensionless.
    CSF_max : float, optional
        Maximum companion star fraction, dimensionless mean companion
        count (not bounded by 1). Default 3.
    q_power : float, optional
        Mass-ratio power-law index, dimensionless. Default -0.4.
    q_min : float, optional
        Minimum mass ratio m_comp/m_prim, dimensionless. Default 0.01.
    companion_max : bool, optional
        If True, cap companion counts at CSF_max. Default False.
    binary_only_mass_max : float, optional
        Primary mass in solar masses (Msun) at and below which systems
        are binaries only. Default 0.08 Msun.
    """
    def __init__(self, MF_A, MF_B, MF_M0, MF_k,
                 CSF_A, CSF_B, CSF_M0, CSF_k,
                 CSF_max=3, q_power=-0.4, q_min=0.01, companion_max=False,
                 binary_only_mass_max=H_BURNING_MASS):
        super(MultiplicityLogistic, self).__init__(
            MF_amp=1.0, MF_power=0.0, CSF_amp=1.0, CSF_power=0.0,
            CSF_max=CSF_max, q_power=q_power, q_min=q_min,
            companion_max=companion_max,
            binary_only_mass_max=binary_only_mass_max)
        self.MF_A = float(MF_A)
        self.MF_B = float(MF_B)
        self.MF_M0 = float(MF_M0)
        self.MF_k = float(MF_k)
        self.CSF_A = float(CSF_A)
        self.CSF_B = float(CSF_B)
        self.CSF_M0 = float(CSF_M0)
        self.CSF_k = float(CSF_k)

    def multiplicity_fraction(self, mass):
        """
        Multiplicity fraction as a logistic in log primary mass.
        Clipped to [0, 1].

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        mf : float or ndarray
            Multiplicity fraction, dimensionless, in [0, 1].
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        return _logistic_in_logm(
            mass, self.MF_A, self.MF_B, self.MF_M0, self.MF_k,
            clip_min=0.0, clip_max=1.0)

    def companion_star_fraction(self, mass):
        """
        Companion star fraction as a logistic in log primary mass.

        Clipped to [0, CSF_max], raised to at least MF, and set equal
        to MF for primaries at or below ``binary_only_mass_max``.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        csf : float or ndarray
            Companion star fraction, dimensionless mean companion
            count (not bounded by 1). Python float if ``mass`` is
            scalar, ndarray otherwise.
        """
        return_scalar = np.isscalar(mass)
        mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
        csf = _logistic_in_logm(
            mass_arr, self.CSF_A, self.CSF_B, self.CSF_M0, self.CSF_k,
            clip_min=0.0, clip_max=self.CSF_max)
        mf = _logistic_in_logm(
            mass_arr, self.MF_A, self.MF_B, self.MF_M0, self.MF_k,
            clip_min=0.0, clip_max=1.0)
        csf = np.maximum(csf, mf)
        bd = mass_arr <= self.binary_only_mass_max
        csf[bd] = mf[bd]
        if return_scalar:
            return float(csf[0])
        return csf


class MultiplicityUnresolvedOffner2023(MultiplicityLogistic):
    """
    Unresolved multiplicity from Offner et al. 2023 Table 1, including
    brown dwarfs.

    Citation: Offner, S. S. R., Moe, M., Kratter, K. M., Sadavoy, S. I.,
    Jensen, E. L. N., & Tobin, J. J. 2023, in Protostars and Planets VII,
    ASP Conf. Ser. 534, 275 (`arXiv:2203.10066
    <https://arxiv.org/abs/2203.10066>`_; ADS
    `2023ASPC..534..275O
    <https://ui.adsabs.harvard.edu/abs/2023ASPC..534..275O>`_).
    Table 1 data: Zenodo `10.5281/zenodo.6628915
    <https://doi.org/10.5281/zenodo.6628915>`_.

    The multiplicity fraction and companion frequency are a
    **4-parameter logistic in log-mass** fitted with equal weight to
    the geom-mean MF/CF columns of Table 1::

        f(M) = A + (B - A) / (1 + (M / M0)**(-k))

    with (A, B, M0, k) = (0.14, 0.99, 1.41, 1.25) for MF and
    (0.12, 2.35, 3.57, 0.96) for CSF/CF. The curve is C-infinity
    smooth (not a broken power law), saturates at B ~ 1 for MF so
    A/B stars stay near the Raghavan/MDS/Sana points, and has a
    low-mass floor A ~ 0.14. Fontanive et al. (2018) 8 ± 6% sits
    ~0.07 below the curve (~15%), which is consistent with the
    Burgasser/Close BD points and within ~1–2σ of Fontanive. MF is
    clipped to [0, 1]. Below 0.08 Msun, CSF = MF (binaries only;
    THF is tiny).

    Companion assignment vs Table 1
    -------------------------------
    Offner et al. 2023 (text above Table 1): BD primaries include all
    BD companions; FGKM MS statistics include only MS companions with
    M_comp > 0.075 Msun; OBA include MS companions above q > 0.1.
    Table 1 stellar MF/CF therefore exclude BD companions. SPISEA still
    draws companions down to ``q_min`` (default 0.01), so brown-dwarf
    secondaries of stellar primaries are generated. The solar-type
    BD-companion fraction is only ≈ 4% (BD desert at a < 0.5 au), so
    the integrated stellar MF is affected very little. Do not interpret
    the simulated stellar-primary MF as a stellar-companion-only
    statistic.

    Mass-ratio draws use an **error-weighted logistic in log-mass**
    fitted to Table 1 γ_trunc (1–100 au)::

        γ(M) = A + (B - A) / (1 + (M / M0)**(-k))

    with (A, B, M0, k) = (6.6, −1.77, 0.0651, 0.629). BD companions
    are still more equal-mass than solar-type companions. The
    err-weighted fit undershoots Fontanive 4.8 ± 2.2 (~3.3 at
    0.033 Msun).

    Characteristic separation μ(a) is a **smooth broken power law**
    in log10(a) vs log10(M) (s = 0.1 dex, FGK-pulled), C-infinity
    via a stable logcosh. σ(log10 a) is a 2-parameter logistic
    pinned at 0.7 / 1.5. See :meth:`log_a_mean` and
    :meth:`sigma_log_a`. Resolved draws use those as loc / scale
    of a truncated lognormal.

    This class is opt-in; it does not change the SPISEA v2.5
    :class:`MultiplicityUnresolved` default.

    Parameters
    ----------
    CSF_max : float, optional
        Maximum companion star fraction, dimensionless mean companion
        count (not bounded by 1). Default 3.
    q_power : float, optional
        Fallback mass-ratio power-law index, dimensionless. Ignored for
        draws when primary mass is provided (the γ logistic is used);
        used by ``random_q(x)`` with no mass. Default 0.2.
    q_min : float, optional
        Minimum mass ratio m_comp/m_prim, dimensionless, in [q_min, 1].
        Default 0.01.
    companion_max : bool, optional
        If True, cap companion counts at CSF_max. Default False.
    binary_only_mass_max : float, optional
        Primary mass in solar masses (Msun) at and below which systems
        are binaries only (CSF = MF, at most one companion). Default
        0.08 Msun.
    """
    def __init__(self, CSF_max=3, q_power=0.2, q_min=0.01,
                 companion_max=False, binary_only_mass_max=H_BURNING_MASS):
        super(MultiplicityUnresolvedOffner2023, self).__init__(
            MF_A=OFFNER2023_MF_A,
            MF_B=OFFNER2023_MF_B,
            MF_M0=OFFNER2023_MF_M0,
            MF_k=OFFNER2023_MF_K,
            CSF_A=OFFNER2023_CSF_A,
            CSF_B=OFFNER2023_CSF_B,
            CSF_M0=OFFNER2023_CSF_M0,
            CSF_k=OFFNER2023_CSF_K,
            CSF_max=CSF_max, q_power=q_power, q_min=q_min,
            companion_max=companion_max,
            binary_only_mass_max=binary_only_mass_max)
        # Table 1/2 data that was fit; evaluation uses the smooth functions.
        self.q_mass = np.array(OFFNER2023_Q_MASS, dtype=float)
        self.q_gamma = np.array(OFFNER2023_Q_GAMMA, dtype=float)

    def q_power_at_mass(self, mass):
        """
        Mass-ratio power-law index γ(M), P(q) ∝ q^γ.

        Error-weighted logistic in log-mass fitted to Table 1
        γ_trunc. Not an interpolation of the Table 1 knots.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        gamma : float or ndarray
            Mass-ratio power-law index γ, dimensionless.
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        return _logistic_in_logm(
            mass, OFFNER2023_Q_A, OFFNER2023_Q_B,
            OFFNER2023_Q_M0, OFFNER2023_Q_K)

    def log_a_mean(self, mass):
        """
        Characteristic log10(a/AU) from the smooth broken power law.

        FGK-pulled, s = 0.1 dex, C-infinity (stable logcosh).
        Linear-space a is clipped to 0.1 AU.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        log_a_mean : float or ndarray
            Characteristic log10(a / 1 AU) in dex (not ln, not AU).
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        return _smooth_broken_loglog(
            mass, OFFNER2023_A_MUP, OFFNER2023_A_MP,
            OFFNER2023_A_ALPHAL, OFFNER2023_A_ALPHAR, OFFNER2023_A_S,
            a_min=OFFNER2023_A_MIN)

    def a_mean(self, mass):
        """
        Characteristic μ(a) in AU, ``10 ** log_a_mean(mass)``.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        a_mean : float or ndarray
            Characteristic separation μ(a) in AU. Python float if
            ``mass`` is scalar, ndarray otherwise.
        """
        log_a = self.log_a_mean(mass)
        if np.isscalar(log_a):
            return 10.0 ** log_a
        return 10.0 ** np.asarray(log_a, dtype=float)

    def sigma_log_a(self, mass):
        """
        σ(log10 a) from a 2-parameter logistic in log-mass.

        Floors/ceilings pinned at 0.7 / 1.5; clipped to ≥ 0.1.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        sigma_log_a : float or ndarray
            Standard deviation of log10(a / 1 AU), in dex.
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        return _logistic_in_logm(
            mass, OFFNER2023_SIG_A, OFFNER2023_SIG_B,
            OFFNER2023_SIG_M0, OFFNER2023_SIG_K, clip_min=0.1)

    def _q_values_for_primaries(self, prim_subset, n_comp, rng):
        """
        Draw mass-dependent q for every primary (BD and stellar).

        Parameters
        ----------
        prim_subset : array_like
            Primary masses in solar masses (Msun) for this companion-count
            group.
        n_comp : int
            Number of companions per primary, integer count.
        rng : numpy.random.Generator
            Random generator used for inverse-CDF q draws.

        Returns
        -------
        q_values : ndarray
            Companion mass ratios m_comp/m_prim, dimensionless, in
            [q_min, 1]. Shape (len(prim_subset), n_comp).
        """
        return self.random_q(
            rng.random((len(prim_subset), n_comp)), mass=prim_subset)


class MultiplicityResolvedOffner2023(MultiplicityUnresolvedOffner2023,
                                     _ResolvedOrbitalMixin):
    """
    Resolved Offner et al. 2023 multiplicity: Table 1 MF/CF plus
    mass-dependent separations.

    Separations are drawn from a truncated lognormal in log10(a/AU)
    with loc = :meth:`log_a_mean` (smooth broken power law, s = 0.1
    dex, FGK-pulled) and scale = :meth:`sigma_log_a` (2-parameter
    logistic pinned at 0.7 / 1.5). Brown-dwarf binaries peak near a
    few AU. Truncation is 0.01–2000 AU, same as
    :class:`MultiplicityResolvedDK`.

    Eccentricity and Keplerian angles follow Duchêne & Kraus (2013),
    same as :class:`MultiplicityResolvedDK`.

    Opt-in; does not replace :class:`MultiplicityResolvedDK`.

    Parameters
    ----------
    CSF_max : float, optional
        Maximum companion star fraction, dimensionless mean companion
        count (not bounded by 1). Default 3.
    q_power : float, optional
        Fallback mass-ratio power-law index, dimensionless. Ignored for
        draws when primary mass is provided (the γ logistic is used).
        Default 0.2.
    q_min : float, optional
        Minimum mass ratio m_comp/m_prim, dimensionless. Default 0.01.
    companion_max : bool, optional
        If True, cap companion counts at CSF_max. Default False.
    binary_only_mass_max : float, optional
        Primary mass in solar masses (Msun) at and below which systems
        are binaries only. Default 0.08 Msun.
    """
    def __init__(self, **kwargs):
        super(MultiplicityResolvedOffner2023, self).__init__(**kwargs)
        # Table 1/2 data that was fit; draws use log_a_mean / sigma_log_a.
        self.sep_mass = np.array(OFFNER2023_SEP_MASS, dtype=float)
        self.sep_mu_au = np.array(OFFNER2023_SEP_MU_AU, dtype=float)
        self.sep_sig_mass = np.array(OFFNER2023_SEP_SIG_MASS, dtype=float)
        self.sep_sig = np.array(OFFNER2023_SEP_SIG, dtype=float)

    def log_semimajoraxis(self, mass):
        """
        Draw log10(a/AU) from a mass-dependent truncated lognormal.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        log_semimajoraxis : ndarray
            Drawn log10(a / 1 AU) in dex (not ln, not AU), truncated
            so a is between 0.01 AU and 2000 AU.
        """
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        log_a_mean = np.atleast_1d(np.asarray(self.log_a_mean(mass), dtype=float))
        log_a_std = np.atleast_1d(np.asarray(self.sigma_log_a(mass), dtype=float))

        log_a_lower = np.log10(0.01)
        log_a_upper = np.log10(2000)
        a_lower_std = (log_a_lower - log_a_mean) / log_a_std
        a_upper_std = (log_a_upper - log_a_mean) / log_a_std
        return truncnorm.rvs(a_lower_std, a_upper_std,
                             loc=log_a_mean, scale=log_a_std)


# Convenience alias; unresolved Table 1 model is the usual opt-in object.
MultiplicityOffner2023 = MultiplicityUnresolvedOffner2023


def _two_point_powerlaw(mass_1, y_1, mass_2, y_2):
    """
    Amplitude and power for y = A * M**alpha through two (M, y) points.

    Parameters
    ----------
    mass_1, mass_2 : float
        Primary masses in solar masses (Msun) of the two anchor points.
        Must be positive and distinct.
    y_1, y_2 : float
        Ordinate values at ``mass_1`` and ``mass_2``. Units match the
        fitted quantity (dimensionless for MF/γ, mean companion count
        for CSF, AU for characteristic a, dex for σ).

    Returns
    -------
    amp : float
        Power-law amplitude A, in units of y / Msun**alpha.
    power : float
        Power-law index alpha, dimensionless.
    """
    power = np.log(y_2 / y_1) / np.log(mass_2 / mass_1)
    amp = y_1 / (mass_1 ** power)
    return amp, power


def _piecewise_powerlaw(mass, mass_limits, amps, powers, clip_min=None,
                        clip_max=None):
    """
    Evaluate y = A_i * M**alpha_i on mass segments.

    Segment i applies for mass_limits[i] <= M < mass_limits[i+1].
    The first segment also covers M below the lowest limit; the last
    segment is closed on the right.

    Parameters
    ----------
    mass : float or array_like
        Primary mass in solar masses (Msun).
    mass_limits : array_like
        Segment edges in solar masses (Msun), length N+1, strictly
        increasing.
    amps : array_like
        Length-N amplitudes A_i, in units of y / Msun**alpha_i.
    powers : array_like
        Length-N power-law indices alpha_i, dimensionless.
    clip_min, clip_max : float or None, optional
        Optional lower/upper clips on y, in the same units as y.
        ``None`` means no clip on that side.

    Returns
    -------
    y : float or ndarray
        Piecewise power-law value, in the same units as
        ``amps * mass**powers``. Python float if ``mass`` is scalar,
        ndarray otherwise.
    """
    return_scalar = np.isscalar(mass)
    mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
    out = np.empty(mass_arr.shape, dtype=float)
    nseg = len(amps)
    for i in range(nseg):
        lo = mass_limits[i]
        hi = mass_limits[i + 1]
        if i == 0:
            mask = mass_arr < hi
        elif i == nseg - 1:
            mask = mass_arr >= lo
        else:
            mask = (mass_arr >= lo) & (mass_arr < hi)
        out[mask] = amps[i] * np.power(mass_arr[mask], powers[i])
    if clip_min is not None:
        out = np.maximum(out, clip_min)
    if clip_max is not None:
        out = np.minimum(out, clip_max)
    if return_scalar:
        return float(out[0])
    return out


def _logistic_in_logm(mass, A, B, M0, k, clip_min=None, clip_max=None):
    """
    Evaluate y = A + (B - A) / (1 + (M / M0)**(-k)).

    This is a logistic in log-mass: as M -> 0, y -> A; as M -> inf,
    y -> B. Masses <= 0 are mapped to the low-mass asymptote A.

    Parameters
    ----------
    mass : float or array_like
        Primary mass in solar masses (Msun).
    A, B : float
        Low-mass and high-mass asymptotes, in the same units as y.
        Dimensionless for MF and γ; mean companion count for CSF;
        dex for σ(log10 a).
    M0 : float
        Characteristic mass in solar masses (Msun).
    k : float
        Logistic slope, dimensionless.
    clip_min, clip_max : float or None, optional
        Optional lower/upper clips on y, in the same units as y.
        ``None`` means no clip on that side.

    Returns
    -------
    y : float or ndarray
        Logistic value in the same units as ``A`` and ``B``.
        Python float if ``mass`` is scalar, ndarray otherwise.
    """
    return_scalar = np.isscalar(mass)
    mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
    out = np.full(mass_arr.shape, float(A), dtype=float)
    positive = mass_arr > 0
    if np.any(positive):
        m = mass_arr[positive]
        out[positive] = A + (B - A) / (1.0 + np.power(m / M0, -k))
    if clip_min is not None:
        out = np.maximum(out, clip_min)
    if clip_max is not None:
        out = np.minimum(out, clip_max)
    if return_scalar:
        return float(out[0])
    return out


def _logcosh(x):
    """
    Numerically stable log(cosh(x)).

    Uses |x| + log1p(exp(-2|x|)) - log(2) rather than np.log(np.cosh(x)),
    which overflows for |x| ≳ 700.

    Parameters
    ----------
    x : float or array_like
        Argument of cosh, dimensionless (for the smooth broken power
        law this is v/s, with v and s in dex).

    Returns
    -------
    logcosh_x : float or ndarray
        log(cosh(x)), dimensionless. Same shape as ``x``.
    """
    ax = np.abs(np.asarray(x, dtype=float))
    return ax + np.log1p(np.exp(-2.0 * ax)) - np.log(2.0)


def _smooth_broken_loglog(mass, mup, Mp, alpha_L, alpha_R, s, a_min=0.1):
    """
    Smooth broken power law in log10(a) vs log10(M).

        v  = log10(M / Mp)
        yp = log10(mup)
        log10(a) = yp + 0.5*(αL+αR)*v + 0.5*(αR-αL)*s * logcosh(v/s)

    ``s`` is the smoothing scale in dex (C-infinity; logcosh). Masses
    <= 0 are mapped to 1e-8 Msun. The linear-space value is clipped
    to ``a_min``.

    Parameters
    ----------
    mass : float or array_like
        Primary mass in solar masses (Msun).
    mup : float
        Characteristic separation at the break mass, in AU.
    Mp : float
        Break mass in solar masses (Msun).
    alpha_L, alpha_R : float
        Power-law indices below and above ``Mp``, dimensionless.
    s : float
        Smoothing scale in dex of log10(M / 1 Msun).
    a_min : float, optional
        Minimum linear-space separation in AU. Default 0.1 AU.

    Returns
    -------
    log_a : float or ndarray
        log10(a / 1 AU) in dex (not ln, not AU). Python float if
        ``mass`` is scalar, ndarray otherwise.
    """
    return_scalar = np.isscalar(mass)
    mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
    m = np.where(mass_arr > 0, mass_arr, 1e-8)
    v = np.log10(m / float(Mp))
    yp = np.log10(float(mup))
    log_a = (yp
             + 0.5 * (alpha_L + alpha_R) * v
             + 0.5 * (alpha_R - alpha_L) * s * _logcosh(v / s))
    a = np.maximum(10.0 ** log_a, float(a_min))
    log_a = np.log10(a)
    if return_scalar:
        return float(log_a[0])
    return log_a


def _q_from_powerlaw(x, q_pow, q_min):
    """
    Inverse CDF of P(q) ∝ q**q_pow for q_min <= q <= 1.

    ``q_pow`` may be a scalar or an array broadcastable to ``x``.
    The q_pow = -1 (b = 0) limit is q = q_min**(1 - x).

    Parameters
    ----------
    x : float or array_like
        Uniform random draw, dimensionless, in [0, 1].
    q_pow : float or array_like
        Mass-ratio power-law index γ, dimensionless. Broadcastable
        to ``x``.
    q_min : float
        Minimum mass ratio m_comp/m_prim, dimensionless, in (0, 1].

    Returns
    -------
    q : ndarray
        Companion mass ratio m_comp/m_prim, dimensionless, in
        [q_min, 1]. Same shape as the broadcast of ``x`` and ``q_pow``.
    """
    x = np.asarray(x, dtype=float)
    q_pow = np.asarray(q_pow, dtype=float)
    if x.ndim > q_pow.ndim:
        q_pow = q_pow.reshape(q_pow.shape + (1,) * (x.ndim - q_pow.ndim))
    b = 1.0 + q_pow
    b, x = np.broadcast_arrays(b, x)
    out = np.empty(x.shape, dtype=float)
    near_zero = np.abs(b) < 1e-12
    ok = ~near_zero
    if np.any(near_zero):
        out[near_zero] = q_min ** (1.0 - x[near_zero])
    if np.any(ok):
        out[ok] = (x[ok] * (1.0 - q_min ** b[ok]) + q_min ** b[ok]) ** (1.0 / b[ok])
    return out


def _offner2023_table1_geom_mass(m_lo, m_hi):
    """
    Geometric-mean primary mass of a Table 1 M1 interval.

    Parameters
    ----------
    m_lo, m_hi : float
        Low and high edges of the Table 1 primary-mass interval, in
        solar masses (Msun). Must be positive.

    Returns
    -------
    mass : float
        Geometric mean sqrt(m_lo * m_hi) in solar masses (Msun).
    """
    return float(np.sqrt(m_lo * m_hi))

# Table 1/2 data that was fit; the model is the smooth functions above,
# not interpolation of these arrays.
OFFNER2023_Q_MASS = np.array([
    _offner2023_table1_geom_mass(0.019, 0.058),  # Fontanive+2018: 4.8
    0.065,                                       # L/early-T: 2-3 (text)
    _offner2023_table1_geom_mass(0.080, 0.095),  # Close+2003: 3.3
    _offner2023_table1_geom_mass(0.06, 0.15),    # Allen+2007: 1.7
    _offner2023_table1_geom_mass(0.15, 0.30),    # Winters mid-M: 0.7
    _offner2023_table1_geom_mass(0.3, 0.6),      # Winters early-M: 0.1
    _offner2023_table1_geom_mass(0.75, 1.25),    # Raghavan FGK: 0.2
    _offner2023_table1_geom_mass(1.6, 2.4),      # De Rosa A: -1.3
    _offner2023_table1_geom_mass(3.0, 5.0),      # MDS 3-5: -1.0
    _offner2023_table1_geom_mass(5.0, 8.0),      # MDS 5-8: -1.7
    _offner2023_table1_geom_mass(8.0, 17.0),     # MDS 8-17: -1.6
    _offner2023_table1_geom_mass(17.0, 50.0),    # Sana O: -1.4
])
OFFNER2023_Q_GAMMA = np.array([
    4.8, 2.5, 3.3, 1.7, 0.7, 0.1, 0.2, -1.3, -1.0, -1.7, -1.6, -1.4
])

# Table 1 ã_all (au) and Table 2 lognormal μ (au) vs geom-mean M1.
# Table 2 μ is listed where both exist (late-M, early-M, FGK).
# Fit data only; evaluation uses the smooth broken power law.
OFFNER2023_SEP_MASS = np.array([
    _offner2023_table1_geom_mass(0.019, 0.058),  # Fontanive ã_all=2.9
    _offner2023_table1_geom_mass(0.080, 0.095),  # Close ã_all=3.7
    _offner2023_table1_geom_mass(0.075, 0.15),   # Table 2 late-M μ=4
    _offner2023_table1_geom_mass(0.15, 0.30),    # Winters mid-M ã_all=10
    _offner2023_table1_geom_mass(0.3, 0.6),      # Table 2 early-M μ=25
    _offner2023_table1_geom_mass(0.75, 1.25),    # Table 2 FGK μ=40
    _offner2023_table1_geom_mass(1.6, 2.4),      # Moe & Kratter ã_all=32
    _offner2023_table1_geom_mass(3.0, 5.0),      # MDS ã_all=28
    _offner2023_table1_geom_mass(5.0, 8.0),      # MDS ã_all=25
    _offner2023_table1_geom_mass(8.0, 17.0),     # MDS ã_all=23
    _offner2023_table1_geom_mass(17.0, 50.0),    # Sana ã_all=19
])
OFFNER2023_SEP_MU_AU = np.array([
    2.9, 3.7, 4.0, 10.0, 25.0, 40.0, 32.0, 28.0, 25.0, 23.0, 19.0
])
# Table 2 σ_log a at the three published bins. Fit data only;
# evaluation uses the 2-parameter logistic.
OFFNER2023_SEP_SIG_MASS = np.array([
    _offner2023_table1_geom_mass(0.075, 0.15),
    _offner2023_table1_geom_mass(0.3, 0.6),
    _offner2023_table1_geom_mass(0.75, 1.25),
])
OFFNER2023_SEP_SIG = np.array([0.7, 1.3, 1.5])
