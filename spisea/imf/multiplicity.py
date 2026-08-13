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


def _two_point_powerlaw(mass_1, y_1, mass_2, y_2):
    """
    Amplitude and power for y = A * M**alpha through two (M, y) points.
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
    segment is closed on the right. Returns a Python float for scalar
    mass and an ndarray otherwise.
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


def _q_from_powerlaw(x, q_pow, q_min):
    """
    Inverse CDF of P(q) ∝ q**q_pow for q_min <= q <= 1.

    ``q_pow`` may be a scalar or an array broadcastable to ``x``.
    The q_pow = -1 (b = 0) limit is q = q_min**(1 - x).
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
        Primary mass (Msun) at and below which systems are restricted
        to at most one companion (CSF = MF). Default is 0.08.
    
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
        mass : float or numpy array
            Mass of primary star.

        Returns
        -------
        mf : float or numpy array
            Multiplicity Fraction, the fraction of stars at this mass
            that will have one or more companions.
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
        mass : float or numpy array
            Mass of primary star

        Returns
        -------
        csf : float or numpy array
            Companion Star Fraction, the expected number of companions
            for a star at this mass.
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
            Random number between 0 and 1.

        mass : float or array_like, optional
            Primary mass(es). If given, the power-law index is
            ``q_power_at_mass(mass)`` (brown-dwarf vs stellar for
            Lu et al. 2013; mass-dependent for Offner et al. 2023).
            If omitted, ``self.q_pow`` is used for all companions.

        Returns
        -------
        q : float or array_like
            companion mass ratio(s)
        """
        if mass is None:
            return _q_from_powerlaw(x, self.q_pow, self.q_min)
        return _q_from_powerlaw(x, self.q_power_at_mass(mass), self.q_min)

    def random_is_multiple(self, x, MF):
        """
        Helper function: determine if star is in multiple system.
        """
        return x < MF

    def random_companion_count(self, x, CSF, MF):
        """
        Helper function: calculate number of companions.
        """
        # bd stipulation since mf=0
        if MF <= 0:
            return 0

        n_comp = 1 + np.random.poisson((CSF / MF) - 1)
        
        if self.companion_max == True:
            if n_comp > self.CSF_max:
                n_comp = self.CSF_max
            
        return n_comp

    def draw_n_companions(self, mass, CSF, MF, rng):
        """
        Vectorized companion counts for primaries that are already
        identified as multiple.

        Parameters
        ----------
        mass, CSF, MF : array_like
            Primary masses and corresponding CSF and MF values.
        rng : numpy.random.Generator
            Random generator used for the Poisson draw.

        Returns
        -------
        n_comp : ndarray of int
            Number of companions for each multiple primary. Brown-dwarf
            primaries are limited to one companion.
        """
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        CSF = np.atleast_1d(np.asarray(CSF, dtype=float))
        MF = np.atleast_1d(np.asarray(MF, dtype=float))
        n_comp = 1 + rng.poisson((CSF / MF) - 1)
        if self.companion_max:
            n_comp = np.minimum(n_comp, self.CSF_max)
        bd = mass <= self.binary_only_mass_max
        n_comp[bd & (n_comp > 1)] = 1
        return n_comp

    def _q_values_for_primaries(self, prim_subset, n_comp, rng):
        """
        Draw mass ratios for ``n_comp`` companions of each primary.

        The stellar / brown-dwarf split (two separate RNG draws) preserves
        the historical Lu et al. (2013) random sequence used by
        ``imf.calc_multi``.
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
            Primary masses (Msun).
        is_multiple : array_like of bool
            True for primaries drawn as multiple systems.
        CSF, MF : array_like
            Companion star fraction and multiplicity fraction at each
            primary mass.
        rng : numpy.random.Generator
            Random generator.
        mass_min : float
            Minimum companion mass; lighter companions are masked.

        Returns
        -------
        comp_masses : numpy.ma.MaskedArray
            Companion masses, shape (n_primaries, max_n_comp).
        system_masses : ndarray
            Primary plus unmasked companion mass.
        is_multiple : ndarray of bool
            Updated multiplicity flags after masking sub-minimum companions.
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


def _offner2023_table1_geom_mass(m_lo, m_hi):
    """Geometric-mean primary mass of a Table 1 M1 interval."""
    return float(np.sqrt(m_lo * m_hi))


def offner2023_default_segments():
    """
    Piecewise MF/CF power-law segments fitted to Offner et al. 2023 Table 1.

    Each segment is the two-point power law through the geometric-mean
    primary masses of the listed surveys. MF is clipped to [0, 1] at
    evaluation time. Below the hydrogen-burning limit, CSF is forced
    equal to MF (binaries only) regardless of the CF segment.

    Returns
    -------
    mass_limits : ndarray
        Segment edges (Msun), length N+1.
    MF_amps, MF_powers, CSF_amps, CSF_powers : ndarray
        Length-N piecewise coefficients for MF(M) = A * M**alpha
        and CSF(M) = A * M**alpha.
    """
    # Geometric-mean M1 of Table 1 mass bins (Offner et al. 2023, Table 1).
    M_font = _offner2023_table1_geom_mass(0.019, 0.058)   # Fontanive+2018
    M_burg = _offner2023_table1_geom_mass(0.05, 0.08)     # Burgasser 2007
    M_close = _offner2023_table1_geom_mass(0.080, 0.095)  # Close+2003
    M_wlate = _offner2023_table1_geom_mass(0.075, 0.15)   # Winters+2019 late-M
    M_wmid = _offner2023_table1_geom_mass(0.15, 0.30)     # Winters+2019 mid-M
    M_wearly = _offner2023_table1_geom_mass(0.3, 0.6)     # Winters+2019 early-M
    M_ragh = _offner2023_table1_geom_mass(0.75, 1.25)     # Raghavan+2010 FGK
    M_moe = _offner2023_table1_geom_mass(1.6, 2.4)        # Moe & Kratter 2021
    M_mds35 = _offner2023_table1_geom_mass(3.0, 5.0)      # Moe & Di Stefano 2017
    M_mds58 = _offner2023_table1_geom_mass(5.0, 8.0)
    M_mds817 = _offner2023_table1_geom_mass(8.0, 17.0)
    M_sana = _offner2023_table1_geom_mass(17.0, 50.0)     # Sana et al.

    # Published Table 1 MF (fraction) and CF. MF values that were given
    # as percents in the table are converted here.
    mf_font, cf_font = 0.08, 0.08
    mf_burg, cf_burg = 0.15, 0.16
    mf_close, cf_close = 0.19, 0.19
    mf_wlate, cf_wlate = 0.19, 0.21
    mf_wmid, cf_wmid = 0.23, 0.27
    mf_wearly, cf_wearly = 0.30, 0.38
    mf_ragh, cf_ragh = 0.46, 0.60
    mf_moe, cf_moe = 0.68, 0.99
    mf_mds35, cf_mds35 = 0.81, 1.28
    mf_mds58, cf_mds58 = 0.89, 1.55
    mf_mds817, cf_mds817 = 0.93, 1.80
    mf_sana, cf_sana = 0.96, 2.10

    # Mass edges follow BD / late-M / mid-M / FGK / A / B / O regimes
    # in Offner et al. 2023 §2.2. Matching Table 1 is preferred over
    # continuity at the breaks.
    mass_limits = np.array([0.01, 0.05, 0.08, 0.15, 0.30, 1.5, 5.0, 17.0, 150.0])

    mf_pairs = [
        (M_font, mf_font, M_burg, mf_burg),       # late-T/Y BDs
        (M_burg, mf_burg, M_close, mf_close),     # L/early-T / upper BD
        (M_wlate, mf_wlate, M_wmid, mf_wmid),     # late-M
        (M_wmid, mf_wmid, M_wearly, mf_wearly),   # mid-M
        (M_wearly, mf_wearly, M_ragh, mf_ragh),   # early-M + FGK
        (M_moe, mf_moe, M_mds35, mf_mds35),       # A / late-B
        (M_mds58, mf_mds58, M_mds817, mf_mds817), # B
        (M_mds817, mf_mds817, M_sana, mf_sana),   # O
    ]
    cf_pairs = [
        (M_font, cf_font, M_burg, cf_burg),
        (M_burg, cf_burg, M_close, cf_close),
        (M_wlate, cf_wlate, M_wmid, cf_wmid),
        (M_wmid, cf_wmid, M_wearly, cf_wearly),
        (M_wearly, cf_wearly, M_ragh, cf_ragh),
        (M_moe, cf_moe, M_mds35, cf_mds35),
        (M_mds58, cf_mds58, M_mds817, cf_mds817),
        (M_mds817, cf_mds817, M_sana, cf_sana),
    ]

    MF_amps = np.empty(len(mf_pairs))
    MF_powers = np.empty(len(mf_pairs))
    CSF_amps = np.empty(len(cf_pairs))
    CSF_powers = np.empty(len(cf_pairs))
    for i, pair in enumerate(mf_pairs):
        MF_amps[i], MF_powers[i] = _two_point_powerlaw(*pair)
    for i, pair in enumerate(cf_pairs):
        CSF_amps[i], CSF_powers[i] = _two_point_powerlaw(*pair)
    return mass_limits, MF_amps, MF_powers, CSF_amps, CSF_powers


# Table 1 gamma_trunc (1-100 au when available) plus the paper's
# L/early-T value gamma = 2-3 (midpoint 2.5). Masses are geometric
# means of the corresponding Table 1 M1 bins.
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
# Table 2 μ is used where both exist (late-M, early-M, FGK).
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
# Table 2 σ_log a at the three published bins; held constant outside.
OFFNER2023_SEP_SIG_MASS = np.array([
    _offner2023_table1_geom_mass(0.075, 0.15),
    _offner2023_table1_geom_mass(0.3, 0.6),
    _offner2023_table1_geom_mass(0.75, 1.25),
])
OFFNER2023_SEP_SIG = np.array([0.7, 1.3, 1.5])


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
        Segment edges in solar masses, length N+1, strictly increasing.
    MF_amps, MF_powers : array_like
        Length-N amplitudes and powers for the multiplicity fraction.
    CSF_amps, CSF_powers : array_like
        Length-N amplitudes and powers for the companion star fraction.
    CSF_max, q_power, q_min, companion_max, binary_only_mass_max
        Passed to :class:`MultiplicityUnresolved`.
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
        """
        return _piecewise_powerlaw(
            mass, self.mass_limits, self.MF_amps, self.MF_powers,
            clip_min=0.0, clip_max=1.0)

    def companion_star_fraction(self, mass):
        """
        Companion star fraction as a piecewise power law in primary mass.

        Clipped to [0, CSF_max], raised to at least MF, and set equal
        to MF for primaries at or below ``binary_only_mass_max``.
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


class MultiplicityUnresolvedOffner2023(MultiplicityPiecewisePowerLaw):
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

    The multiplicity fraction and companion frequency are piecewise
    power laws fitted to the bias-corrected Table 1 MF and CF at the
    geometric-mean primary mass of each survey bin. Table 1 itself
    does **not** publish MF(M) ∝ M^α; the only tabulated power law is
    γ_trunc for f_q ∝ q^γ on 0.4 < q < 1. Matching Table 1 is preferred
    over continuity at segment boundaries. MF is clipped to [0, 1].
    Below 0.08 Msun, CSF = MF (binaries only; THF is tiny).

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

    Mass-ratio draws use Table 1 γ_trunc (1–100 au when listed),
    interpolated in log mass, so BD companions are much more
    equal-mass than solar-type companions.

    This class is opt-in; it does not change the Lu et al. (2013)
    :class:`MultiplicityUnresolved` default.

    Parameters
    ----------
    CSF_max, q_min, companion_max, binary_only_mass_max
        See :class:`MultiplicityUnresolved`. ``q_power`` is ignored for
        draws when primary mass is provided (Table 1 γ_trunc is used);
        it remains the fallback for ``random_q(x)`` with no mass.
    """
    def __init__(self, CSF_max=3, q_power=0.2, q_min=0.01,
                 companion_max=False, binary_only_mass_max=H_BURNING_MASS):
        mass_limits, MF_amps, MF_powers, CSF_amps, CSF_powers = \
            offner2023_default_segments()
        super(MultiplicityUnresolvedOffner2023, self).__init__(
            mass_limits, MF_amps, MF_powers, CSF_amps, CSF_powers,
            CSF_max=CSF_max, q_power=q_power, q_min=q_min,
            companion_max=companion_max,
            binary_only_mass_max=binary_only_mass_max)
        self.q_mass = np.array(OFFNER2023_Q_MASS, dtype=float)
        self.q_gamma = np.array(OFFNER2023_Q_GAMMA, dtype=float)

    def q_power_at_mass(self, mass):
        """
        Table 1 γ_trunc interpolated in log primary mass.

        Late-T/Y BDs have γ ≈ 4.8 (Fontanive+2018); L/early-T have
        γ ≈ 2–3; solar-type γ_trunc ≈ 0.2 at 1–100 au; massive stars
        have negative γ (low-q companions).
        """
        mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
        mass_clip = np.clip(mass_arr, self.q_mass[0], self.q_mass[-1])
        q_pow = np.interp(np.log10(mass_clip),
                          np.log10(self.q_mass), self.q_gamma)
        if np.isscalar(mass):
            return float(q_pow[0])
        return q_pow

    def _q_values_for_primaries(self, prim_subset, n_comp, rng):
        """Mass-dependent q for every primary (BD and stellar)."""
        return self.random_q(
            rng.random((len(prim_subset), n_comp)), mass=prim_subset)


class MultiplicityResolvedOffner2023(MultiplicityUnresolvedOffner2023,
                                     _ResolvedOrbitalMixin):
    """
    Resolved Offner et al. 2023 multiplicity: Table 1 MF/CF plus
    mass-dependent separations.

    Separations are drawn from a truncated lognormal in log10(a/AU).
    The mean μ(a) interpolates Table 1 ã_all and Table 2 lognormal μ
    (late-M μ=4 au, σ=0.7; early-M μ=25 au, σ=1.3; FGK μ=40 au,
    σ=1.5). Brown-dwarf / late-M binaries peak near 3–4 au (vast
    majority at 1–20 au). σ_log a uses Table 2 and is held at the
    late-M value (0.7) through the BD regime, as Offner et al. 2023
    describe the BD separation distribution as similar to late-M.

    Eccentricity and Keplerian angles follow Duchêne & Kraus (2013),
    same as :class:`MultiplicityResolvedDK`.

    Opt-in; does not replace :class:`MultiplicityResolvedDK`.
    """
    def __init__(self, **kwargs):
        super(MultiplicityResolvedOffner2023, self).__init__(**kwargs)
        self.sep_mass = np.array(OFFNER2023_SEP_MASS, dtype=float)
        self.sep_mu_au = np.array(OFFNER2023_SEP_MU_AU, dtype=float)
        self.sep_sig_mass = np.array(OFFNER2023_SEP_SIG_MASS, dtype=float)
        self.sep_sig = np.array(OFFNER2023_SEP_SIG, dtype=float)

    def log_semimajoraxis(self, mass):
        """
        Draw log10(a/AU) from a mass-dependent truncated lognormal.

        Parameters
        ----------
        mass : array-like
            Primary mass (Msun).

        Returns
        -------
        log_semimajoraxis : ndarray
            log10 of the semimajor axis in AU, truncated to 0.01–2000 AU.
        """
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        mass_clip = np.clip(mass, self.sep_mass[0], self.sep_mass[-1])
        logm = np.log10(mass_clip)
        log_a_mean = np.interp(logm, np.log10(self.sep_mass),
                               np.log10(self.sep_mu_au))
        log_a_std = np.interp(logm, np.log10(self.sep_sig_mass), self.sep_sig)
        log_a_std = np.clip(log_a_std, 0.1, None)

        log_a_lower = np.log10(0.01)
        log_a_upper = np.log10(2000)
        a_lower_std = (log_a_lower - log_a_mean) / log_a_std
        a_upper_std = (log_a_upper - log_a_mean) / log_a_std
        return truncnorm.rvs(a_lower_std, a_upper_std,
                             loc=log_a_mean, scale=log_a_std)


# Convenience alias; unresolved Table 1 model is the usual opt-in object.
MultiplicityOffner2023 = MultiplicityUnresolvedOffner2023
