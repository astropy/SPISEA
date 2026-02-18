import numpy as np
import astropy.modeling
from random import choice

defaultMF_amp = 0.44
defaultMF_power = 0.51
defaultCSF_amp = 0.50
defaultCSF_power = 0.45
defaultCSF_max = 3
defaultq_power = -0.4
defaultq_min = 0.01
default_aMean = 100.0 # log (AU)
default_aSigma = 0.1  # log (AU)

# Eventually we should add in separation properties. (a_mean, a_sigma)

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

    **Companion Star Fraction** -- the expected number of companions in
    a multiple system. The companion star fraction (CSF) also 
    changes with mass and this dependency can be described as
    a power-law::
                
                CSF(mass) = CSF_amp * (mass ** CSF_power)

    The companion star fraction is clipped to some maximum
    value, CSF_max. The actual number of companions is drawn 
    from a Poisson distribution with an expectation value of CSF.

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
    
    """
    def __init__(self, 
                 MF_amp=0.44, MF_power=0.51,
                 CSF_amp=0.50, CSF_power=0.45, CSF_max=3,
                 q_power=-0.4, q_min=0.01, companion_max = False):
         
        self.MF_amp = MF_amp
        self.MF_pow = MF_power
        self.CSF_amp = CSF_amp
        self.CSF_pow = CSF_power
        self.CSF_max = CSF_max
        self.q_pow = q_power
        self.q_min = q_min
        self.companion_max = companion_max

    def multiplicity_fraction(self, mass):
        """
        Given a star's mass, determine the probability that the star is in a
        multiple system (multiplicity fraction = MF).

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
        else:
            mf[mf > 1] = 1

        return mf

    def companion_star_fraction(self, mass):
        """
        Given a star's mass, determine the average number of
        companion stars (companion star fraction = CSF).
        For only binaries (no triples), we will return
        the Multiplicity Fraction back again and set the companion_max = 1.

        Parameters
        ----------
        mass : float or numpy array
            Mass of primary star.

        Returns
        -------
        csf : float or numpy array
            Companion Star Fraction, the expected number of companions
            for a star at this mass.
        """
        return self.multiplicity_fraction(mass)

    def random_q(self, x):
        """
        Generative function for companion mass ratio, equivalent
        to the inverse of the CDF.

            `q = m_compnaion / m_primary`
            `P(q) = q ** beta`    for q_min <= q <= 1

        Parameters
        ----------
        x : float or array_like
            Random number between 0 and 1.

        Returns
        -------
        q : float or array_like
            companion mass ratio(s)
        """
        b = 1.0 + self.q_pow
        q = (x * (1.0 - self.q_min ** b) + self.q_min ** b) ** (1.0 / b)

        return  q

    def get_resolved_companions(self, mass1)
        """
        Function that generates companion masses and orbital
        parameters.

        Not defined for the unresolved class, will be defined in
        resolved subclasses.
        """
        raise NotImplementedError("Function get_resolved_companions is not"
            " defined for MultiplicityUnresolved, only its resolved subclasses.")


    def random_is_multiple(self, x, MF):
        """
        Helper function: determine if star is in multiple system.
        """
        return x < MF

    def random_companion_count(self, x, CSF, MF):
        """
        Helper function: calculate number of companions.
        """
        n_comp = 1 + np.random.poisson((CSF / MF) - 1)
        
        if self.companion_max == True:
            if n_comp > self.CSF_max:
                n_comp = self.CSF_max
            
        return n_comp
    
class MultiplicityResolvedDK(MultiplicityUnresolved):
    """
    Sub-class of MultiplicityUnresolved that adds semimajor axis and eccentricity information 
    for multiple objects from distributions described in Duchene and Kraus 2013
    
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
        
        Parameters
        ----------
        mass : float
            Mass of primary star

        Returns
        -------
        log_semimajoraxis : float
            Log of the semimajor axis/separation between the stars in units of AU
        """
        a_mean_func = astropy.modeling.powerlaws.BrokenPowerLaw1D(amplitude=self.a_amp, x_break=self.a_break, alpha_1=self.a_slope1, alpha_2=self.a_slope2)
        log_a_mean = np.log10(a_mean_func(mass)) #mean log(a)
        log_a_std_func = astropy.modeling.models.Linear1D(slope=self.a_std_slope, intercept=self.a_std_intercept)
        log_a_std = log_a_std_func(np.log10(mass)) #sigma_log(a)
        if mass >= 2.9:
            log_a_std = log_a_std_func(np.log10(2.9)) #sigma_log(a)
        if log_a_std < 0.1:
            log_a_std = 0.1
            
        log_semimajoraxis = np.random.normal(log_a_mean, log_a_std)
        while 10**log_semimajoraxis > 2000 or log_semimajoraxis < -2: #AU
            log_semimajoraxis = np.random.normal(log_a_mean, log_a_std)
            
        return log_semimajoraxis
    
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
            companion mass ratio(s)
        """
        e = np.sqrt(x)
        
        return e

    def get_resolved_companions(self, mass1):
        """
        Function that generates companion masses and orbital
        parameters.

        Parameters
        ----------
        mass1 : float or array_like
            initial masses of primary stars for singles + systems

        Returns
        companions_table : astropy Table with companion initial parameters
        """
        # TODO: bring together the number, mass, a, and e generators here
    
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


class Multiplicity_MoeDiStefano(MultiplicityUnresolved):
    def __init__(self, **kwargs):
        """
        Code modified from Max Moe's IDL function and from the python
        implementation of Moe's code in COSMIC.
        https://raw.githubusercontent.com/COSMIC-PopSynth/COSMIC/refs/heads/develop/src/cosmic/sample/sampler/multidim.py

        Parameters
        ----------
        kwargs
        """

        # Setup the allowed ranges for M1 (primary mass) and M2 (companion masses)
        self.M1min = 0.08
        self.M2min = 0.08
        self.M1max = 150.0
        self.M2max = 150.0

        # Setup the allowed ranges for the orbital periods.
        self.porb_lo = 0.15   # lower limit in log(Period [days])
        self.porb_hi = 8.0    # upper limit in log(Period [days])

        # Tabulate probably density functions of periods,
        # mass ratios, and eccentricities based on
        # analytic fits to corrected binary star populations.
        numM1 = 101

        # use binwidths to maintain structure of original array
        # default size is: numlogP=158
        bwlogP = 0.05
        numq = 91
        nume = 100

        # Vector of primary masses M1 (Msun), logarithmic orbital period P (days),
        # mass ratios q = Mcomp/M1, and eccentricities e
        #
        # 0.8 < M1 < 40 (where we have statistics corrected for selection effects)
        M1_lo = 0.8
        M1_hi = 40
        self.M1v = np.logspace(np.log10(M1_lo), np.log10(M1_hi), numM1)

        # 0.15 < log P < 8.0
        # or use user specified values
        log10_porb_lo = self.porb_lo
        log10_porb_hi = self.porb_hi
        self.logPv = np.arange(log10_porb_lo, log10_porb_hi + bwlogP, bwlogP)
        numlogP = len(self.logPv)

        # 0.10 < q < 1.00
        q_lo = 0.1
        q_hi = 1.0
        self.qv = np.linspace(q_lo, q_hi, numq)

        # 0.0001 < e < 0.9901
        # set minimum to non-zero value to avoid numerical errors
        e_lo = 0.0
        e_hi = 0.99
        self.ev = np.linspace(e_lo, e_hi, nume) + 0.0001
        # Note that companions outside this parameter space (e.g., q < 0.1,
        # log P (days) > 8.0 are not constrained in M+D16 and therefore
        # not considered.

        # Distribution functions - define here, but evaluate within for loops.

        # Frequency of companions with q > 0.1 per decade of orbital period.
        # Bottom panel in Fig. 37 of M+D17
        flogP_sq = np.zeros([numlogP, numM1])

        # Given M1 and P, the cumulative distribution of mass ratios q
        self.cumqdist = np.zeros([numq, numlogP, numM1])

        # Given M1 and P, the cumulative distribution of eccentricities e
        self.cumedist = np.zeros([nume, numlogP, numM1])

        # Given M1 and P, the probability that the companion
        # is a member of the inner binary (currently an approximation).
        # 100% for log P < 1.5, decreases with increasing P
        probbin = np.zeros([numlogP, numM1])

        # Given M1, the cumulative period distribution of the inner binary
        # Normalized so that max(cumPbindist) = total binary frac. (NOT unity)
        cumPbindist = np.zeros([numlogP, numM1])
        # Slope alpha of period distribution across intermediate periods
        # 2.7 - DlogP < log P < 2.7 + DlogP, see Section 9.3 and Eqn. 23.
        # Slightly updated from version 1.
        alpha = 0.018
        DlogP = 0.7

        # Heaviside function for twins with 0.95 < q < 1.00
        H = np.zeros(numq)
        ind = np.where(self.qv >= 0.95)
        H[ind] = 1.0
        H = H / idl_tabulate(self.qv, H)  #normalize so that integral is unity

        # Relevant indices with respect to mass ratio
        indlq = np.where(self.qv >= 0.3)
        indsq = np.where(self.qv < 0.3)
        indq0p3 = np.min(indlq)

        # FILL IN THE MULTIDIMENSIONAL DISTRIBUTION FUNCTIONS
        # Loop through primary mass
        for i in range(0, numM1):
            myM1 = self.M1v[i]
            
            # Twin fraction parameters that are dependent on M1 only; section 9.1
            FtwinlogPle1 = 0.3 - 0.15 * np.log10(myM1)  # Eqn. 6
            logPtwin = 8.0 - myM1  # Eqn. 7a
            if myM1 >= 6.5:
                logPtwin = 1.5  # Eqn. 7b
                
            # Frequency of companions with q > 0.3 at different orbital periods
            # and dependent on M1 only; section 9.3 (slightly modified since v1)
            flogPle1 = (
                    0.020 + 0.04 * np.log10(myM1) + 0.07 * (np.log10(myM1)) ** 2.0
            )  # Eqn. 20
            flogPeq2p7 = (
                    0.039 + 0.07 * np.log10(myM1) + 0.01 * (np.log10(myM1)) ** 2.0
            )  # Eqn. 21
            flogPeq5p5 = (
                    0.078 - 0.05 * np.log10(myM1) + 0.04 * (np.log10(myM1)) ** 2.0
            )  # Eqn. 22

            # Loop through orbital period P
            for j in range(0, numlogP):
                mylogP = self.logPv[j]

                # Given M1 and P, set excess twin fraction; section 9.1 and Eqn. 5
                if mylogP <= 1.0:
                    Ftwin = FtwinlogPle1
                else:
                    Ftwin = FtwinlogPle1 * (1.0 - (mylogP - 1.0) / (logPtwin - 1.0))
                if mylogP >= logPtwin:
                    Ftwin = 0.0

                # Power-law slope gamma_largeq for M1 < 1.2 Msun and various P; Eqn. 9
                if mylogP <= 5.0:
                    gl_1p2 = -0.5
                if mylogP > 5.0:
                    gl_1p2 = -0.5 - 0.3 * (mylogP - 5.0)

                # Power-law slope gamma_largeq for M1 = 3.5 Msun and various P; Eqn. 10
                if mylogP <= 1.0:
                    gl_3p5 = -0.5
                if (mylogP > 1.0) and (mylogP <= 4.5):
                    gl_3p5 = -0.5 - 0.2 * (mylogP - 1.0)
                if (mylogP > 4.5) and (mylogP <= 6.5):
                    gl_3p5 = -1.2 - 0.4 * (mylogP - 4.5)
                if mylogP > 6.5:
                    gl_3p5 = -2.0

                # Power-law slope gamma_largeq for M1 > 6 Msun and various P; Eqn. 11
                if mylogP <= 1.0:
                    gl_6 = -0.5
                if (mylogP > 1.0) and (mylogP <= 2.0):
                    gl_6 = -0.5 - 0.9 * (mylogP - 1.0)
                if (mylogP > 2.0) and (mylogP <= 4.0):
                    gl_6 = -1.4 - 0.3 * (mylogP - 2.0)

                if mylogP > 4.0:
                    gl_6 = -2.0

                # Given P, interpolate gamma_largeq w/ respect to M1 at myM1
                if myM1 <= 1.2:
                    gl = gl_1p2
                if (myM1 > 1.2) and (myM1 <= 3.5):
                    gl = np.interp(
                        np.log10(myM1), np.log10([1.2, 3.5]), [gl_1p2, gl_3p5]
                    )
                if (myM1 > 3.5) and (myM1 <= 6.0):
                    gl = np.interp(np.log10(myM1), np.log10([3.5, 6.0]), [gl_3p5, gl_6])
                if myM1 > 6.0:
                    gl = gl_6

                # Power-law slope gamma_smallq for M1 < 1.2 Msun and all P; Eqn. 13
                gs_1p2 = 0.3

                # Power-law slope gamma_smallq for M1 = 3.5 Msun and various P; Eqn. 14
                if mylogP <= 2.5:
                    gs_3p5 = 0.2
                if (mylogP > 2.5) and (mylogP <= 5.5):
                    gs_3p5 = 0.2 - 0.3 * (mylogP - 2.5)
                if mylogP > 5.5:
                    gs_3p5 = -0.7 - 0.2 * (mylogP - 5.5)

                # Power-law slope gamma_smallq for M1 > 6 Msun and various P; Eqn. 15
                if mylogP <= 1.0:
                    gs_6 = 0.1
                if (mylogP > 1.0) and (mylogP <= 3.0):
                    gs_6 = 0.1 - 0.15 * (mylogP - 1.0)
                if (mylogP > 3.0) and (mylogP <= 5.6):
                    gs_6 = -0.2 - 0.50 * (mylogP - 3.0)
                if mylogP > 5.6:
                    gs_6 = -1.5

                # Given P, interpolate gamma_smallq w/ respect to M1 at myM1
                if myM1 <= 1.2:
                    gs = gs_1p2
                if (myM1 > 1.2) and (myM1 <= 3.5):
                    gs = np.interp(
                        np.log10(myM1), np.log10([1.2, 3.5]), [gs_1p2, gs_3p5]
                    )
                if (myM1 > 3.5) and (myM1 <= 6.0):
                    gs = np.interp(np.log10(myM1), np.log10([3.5, 6.0]), [gs_3p5, gs_6])
                if myM1 > 6.0:
                    gs = gs_6

                # Given Ftwin, gamma_smallq, and gamma_largeq at the specified M1 & P,
                # tabulate the cumulative mass ratio distribution across 0.1 < q < 1.0
                fq = self.qv ** gl  # slope across 0.3 < q < 1.0
                fq = fq / idl_tabulate(self.qv[indlq], fq[indlq])  # normalize to 0.3 < q < 1.0
                fq = fq * (1.0 - Ftwin) + H * Ftwin  # add twins
                fq[indsq] = (fq[indq0p3] * (self.qv[indsq] / 0.3) ** gs)  # slope across 0.1 < q < 0.3
                cumfq = np.cumsum(fq) - fq[0]  # cumulative distribution
                cumfq = cumfq / np.max(cumfq)  # normalize cumfq(q=1.0) = 1
                self.cumqdist[:, j, i] = cumfq  # save to grid

                # Given M1 and P, q_factor is the ratio of all binaries 0.1 < q < 1.0
                # to those with 0.3 < q < 1.0
                q_factor = idl_tabulate(self.qv, fq)

                # Given M1 & P, calculate power-law slope eta of eccentricity dist.
                if mylogP >= 0.7:
                    # For log P > 0.7 use fits in Section 9.2.
                    # Power-law slope eta for M1 < 3 Msun and log P > 0.7
                    eta_3 = 0.6 - 0.7 / (mylogP - 0.5)  # Eqn. 17
                    # Power-law slope eta for M1 > 7 Msun and log P > 0.7
                    eta_7 = 0.9 - 0.2 / (mylogP - 0.5)  # Eqn. 18
                else:
                    # For log P < 0.7, set eta to fitted values at log P = 0.7
                    eta_3 = -2.9
                    eta_7 = -0.1

                # Given P, interpolate eta with respect to M1 at myM1
                if myM1 <= 3.0:
                    eta = eta_3
                if (myM1 > 3.0) and (myM1 <= 7.0):
                    eta = np.interp(
                        np.log10(myM1), np.log10([3.0, 7.0]), [eta_3, eta_7]
                    )
                if myM1 > 7.0:
                    eta = eta_7

                # Given eta at the specified M1 and P, tabulate eccentricity distribution
                if 10 ** mylogP <= 2.0:
                    # For P < 2 days, assume all systems are close to circular
                    # For adopted ev (spacing and minimum value), eta = -3.2 satisfies this
                    fe = self.ev ** (-3.2)
                else:
                    fe = self.ev ** eta
                    e_max = 1.0 - (10 ** mylogP / 2.0) ** (
                            -2.0 / 3.0
                    )  # maximum eccentricity for given P
                    ind = np.where(self.ev >= e_max)
                    fe[ind] = 0.0  # set dist. = 0 for e > e_max
                    # Assume e dist. has power-law slope eta for 0.0 < e / e_max < 0.8 and
                    # then linear turnover between 0.8 < e / e_max < 1.0 so that dist.
                    # is continuous at e / e_max = 0.8 and zero at e = e_max
                    ind = np.where((self.ev >= 0.8 * e_max) & (self.ev <= 1.0 * e_max))
                    ind_cont = np.min(ind) - 1
                    fe[ind] = np.interp(
                        self.ev[ind], [0.8 * e_max, 1.0 * e_max], [fe[ind_cont], 0.0]
                    )

                cumfe = np.cumsum(fe) - fe[0]  # cumulative distribution
                cumfe = cumfe / np.max(cumfe)  # normalize cumfe(e=e_max) = 1
                self.cumedist[:, j, i] = cumfe  # save to grid

                # Given constants alpha and DlogP and
                # M1 dependent values flogPle1, flogPeq2p7, and flogPeq5p5,
                # calculate frequency flogP of companions with q > 0.3 per decade
                # of orbital period at given P (Section 9.3 and Eqn. 23)
                if mylogP <= 1.0:
                    flogP = flogPle1
                if (mylogP > 1.0) and (mylogP <= 2.7 - DlogP):
                    flogP = flogPle1 + (mylogP - 1.0) / (1.7 - DlogP) * (
                            flogPeq2p7 - flogPle1 - alpha * DlogP
                    )
                if (mylogP > 2.7 - DlogP) and (mylogP <= 2.7 + DlogP):
                    flogP = flogPeq2p7 + alpha * (mylogP - 2.7)
                if (mylogP > 2.7 + DlogP) and (mylogP <= 5.5):
                    flogP = (
                            flogPeq2p7
                            + alpha * DlogP
                            + (mylogP - 2.7 - DlogP)
                            / (2.8 - DlogP)
                            * (flogPeq5p5 - flogPeq2p7 - alpha * DlogP)
                    )
                if mylogP > 5.5:
                    flogP = flogPeq5p5 * np.exp(-0.3 * (mylogP - 5.5))

                # Convert frequency of companions with q > 0.3 to frequency of
                # companions with q > 0.1 according to q_factor; save to grid
                flogP_sq[j, i] = flogP * q_factor

                # Calculate prob. that a companion to M1 with period P is the
                # inner binary.  Currently this is an approximation.
                # 100% for log P < 1.5
                # For log P > 1.5 adopt functional form that reproduces M1 dependent
                # multiplicity statistics in Section 9.4, including a
                # 41% binary star faction (59% single star fraction) for M1 = 1 Msun and
                # 96% binary star fraction (4% single star fraction) for M1 = 28 Msun
                if mylogP <= 1.5:
                    probbin[j, i] = 1.0
                else:
                    probbin[j, i] = (
                            1.0 - 0.11 * (mylogP - 1.5) ** 1.43 * (myM1 / 10.0) ** 0.56
                    )
                if probbin[j, i] <= 0.0:
                    probbin[j, i] = 0.0

            # Given M1, calculate cumulative binary period distribution
            mycumPbindist = (
                    np.cumsum(flogP_sq[:, i] * probbin[:, i])
                    - flogP_sq[0, i] * probbin[0, i]
            )

            # Normalize so that max(cumPbindist) = total binary star fraction (NOT 1)
            mycumPbindist = (
                    mycumPbindist
                    / np.max(mycumPbindist)
                    * idl_tabulate(self.logPv, flogP_sq[:, i] * probbin[:, i])
            )

            self.cumPbindist[:, i] = mycumPbindist  #save to grid

        return

    def multiplicity_fraction(self, mass):
        """
        Return the multiplicity fraction for this mass bin.

        Parameters
        ----------
        mass : float or np.array
            Single or list of primary masses.

        Returns
        -------
            Single or list of multiplicity fractions for input star.
        """
        # Full primary mass vector across 0.08 < M1 < 150

        # Find index of M1v that is closest to myM1.
        #     For M1 = 40 - 150 Msun, adopt binary statistics of M1 = 40 Msun.
        #     or M1 = 0.08 - 0.8 Msun, adopt P and e dist of M1 = 0.8Msun,
        #     scale and interpolate the companion frequencies so that the
        #     binary star fraction of M1 = 0.08 Msun primaries is zero,
        #     and truncate the q distribution so that q > q_min = 0.08/M1
        indM1 = np.where(abs(mass - self.M1v) == min(abs(mass - self.M1v)))
        indM1 = indM1[0]

        # Given M1, determine cumulative binary period distribution
        mycumPbindist_flat = (self.cumPbindist[:, indM1]).flatten()
        # If M1 < 0.8 Msun, rescale to appropriate binary star fraction
        if(mass <= 0.8):
            mycumPbindist_flat = mycumPbindist_flat * np.interp(np.log10(myM1), np.log10([0.08, 0.8]), [0.0, 1.0])

        # Given M1, determine the binary star fraction
        mybinfrac = np.max(mycumPbindist_flat)

        return mybinfrac

    def companion_star_fraction(self, mass):
        """
        For Moe and DiStefano, there are no triples.

        Parameters
        ----------
        mass : float or np.array
            Single or list of primary masses.

        Returns
        -------
            Companion star fraction for each star in the input list.
        """

        return np.minimum(1.0, self.multiplicity_fraction(mass))

    def random_q(self, x, mass1):
        # Full primary mass vector across 0.08 < M1 < 150
        myM1 = np.asarray(mass1)

        # Sadly, we have to do this as a for loop

        # Find index of M1v that is closest to myM1.
        #     For M1 = 40 - 150 Msun, adopt binary statistics of M1 = 40 Msun.
        #     or M1 = 0.08 - 0.8 Msun, adopt P and e dist of M1 = 0.8Msun,
        #     scale and interpolate the companion frequencies so that the
        #     binary star fraction of M1 = 0.08 Msun primaries is zero,
        #     and truncate the q distribution so that q > q_min = 0.08/M1
        indM1 = np.where(abs(myM1 - self.M1v) == min(abs(myM1 - self.M1v)))
        indM1 = indM1[0]

        # Given M1, determine cumulative binary period distribution
        mycumPbindist_flat = (self.cumPbindist[:, indM1]).flatten()
        # If M1 < 0.8 Msun, rescale to appropriate binary star fraction
        if(myM1 <= 0.8):
            mycumPbindist_flat = mycumPbindist_flat * np.interp(np.log10(myM1), np.log10([0.08, 0.8]), [0.0, 1.0])

        # Given M1, determine the binary star fraction
        mybinfrac = np.max(mycumPbindist_flat)

        # Generate random number myrand between 0 and 1
        myrand = x

        # Given myrand, select P and corresponding index in logPv
        mylogP = np.interp(myrand, mycumPbindist_flat, self.logPv)
        indlogP = np.where(abs(mylogP - self.logPv) == min(abs(mylogP - self.logPv)))
        indlogP = indlogP[0]

        # Given M1 & P, determine mass ratio distribution.
        # If M1 < 0.8 Msun, truncate q distribution and consider
        # only mass ratios q > q_min = 0.08 / M1
        mycumqdist = self.cumqdist[:, indlogP, indM1].flatten()
        if (myM1 < 0.8):
            q_min = 0.08 / myM1
            # Calculate cumulative probability at q = q_min
            cum_qmin = np.interp(q_min, self.qv, mycumqdist)
            # Rescale and renormalize cumulative distribution for q > q_min
            mycumqdist = mycumqdist - cum_qmin
            mycumqdist = mycumqdist / max(mycumqdist)
            # Set probability = 0 where q < q_min
            indq = np.where(self.qv <= q_min)
            mycumqdist[indq] = 0.0

        # Given M1 & P, select q from cumulative mass ratio distribution
        myq = np.interp(np.random.rand(), mycumqdist, self.qv)

        return myq

    def random_keplarian_parameter(self, mass1, seed=0):
        """
        Given a primary mass, generate companions and return the
        Keplerian orbital parameters for those companion.
        
        Parameters
        ----------
        mass1 : float
            Primary mass.
        
        seed : int
            Default = 0

        Returns
        -------

        """
        # Step 2
        # Implement Monte Carlo method / random number generator to select
        # single stars and binaries from the grids of distributions

        # Create vector for PRIMARY mass function, which is the mass distribution
        # of single stars and primaries in binaries.
        # This is NOT the IMF, which is the mass distribution of single stars,
        # primaries in binaries, and secondaries in binaries.

        np.random.seed(seed)

        mass_singles = 0.0
        mass_binaries = 0.0
        n_singles = 0
        n_binaries = 0
        primary_mass_list = []
        secondary_mass_list = []
        single_mass_list = []
        porb_list = []
        ecc_list = []
        binfrac_list = []

        # Full primary mass vector across 0.08 < M1 < 150
        myM1 = mass1

        # Find index of M1v that is closest to myM1.
        #     For M1 = 40 - 150 Msun, adopt binary statistics of M1 = 40 Msun.
        #     or M1 = 0.08 - 0.8 Msun, adopt P and e dist of M1 = 0.8Msun,
        #     scale and interpolate the companion frequencies so that the
        #     binary star fraction of M1 = 0.08 Msun primaries is zero,
        #     and truncate the q distribution so that q > q_min = 0.08/M1
        indM1 = np.where(abs(myM1 - self.M1v) == min(abs(myM1 - self.M1v)))
        indM1 = indM1[0]

        # Given M1, determine cumulative binary period distribution
        mycumPbindist_flat = (self.cumPbindist[:, indM1]).flatten()
        # If M1 < 0.8 Msun, rescale to appropriate binary star fraction
        if(myM1 <= 0.8):
            mycumPbindist_flat = mycumPbindist_flat * np.interp(np.log10(myM1), np.log10([0.08, 0.8]), [0.0, 1.0])

        # Given M1, determine the binary star fraction
        mybinfrac = np.max(mycumPbindist_flat)

        # Generate random number myrand between 0 and 1
        myrand = np.random.rand()

        # If random number < binary star fraction, generate a binary
        if(myrand < mybinfrac):
            # Given myrand, select P and corresponding index in logPv
            mylogP = np.interp(myrand, mycumPbindist_flat, self.logPv)
            indlogP = np.where(abs(mylogP - self.logPv) == min(abs(mylogP - self.logPv)))
            indlogP = indlogP[0]

            # Given M1 & P, select e from eccentricity distribution
            mye = np.interp(np.random.rand(), self.cumedist[:, indlogP, indM1].flatten(), self.ev)

            # Given M1 & P, determine mass ratio distribution.
            # If M1 < 0.8 Msun, truncate q distribution and consider
            # only mass ratios q > q_min = 0.08 / M1
                mycumqdist = self.cumqdist[:, indlogP, indM1].flatten()
                if(myM1 < 0.8):
                    q_min = 0.08 / myM1
                    # Calculate cumulative probability at q = q_min
                    cum_qmin = np.interp(q_min, self.qv, mycumqdist)
                    # Rescale and renormalize cumulative distribution for q > q_min
                    mycumqdist = mycumqdist - cum_qmin
                    mycumqdist = mycumqdist / max(mycumqdist)
                    # Set probability = 0 where q < q_min
                    indq = np.where(self.qv <= q_min)
                    mycumqdist[indq] = 0.0

                # Given M1 & P, select q from cumulative mass ratio distribution
                myq = np.interp(np.random.rand(), mycumqdist, self.qv)

                if ((myM1 > self.M1min) and (myq * myM1 > self.M2min) and (myM1 < self.M1max) and
                   (myq * myM1 < self.M2max) and (mylogP < self.porb_hi) and (mylogP > self.porb_lo)):
                    primary_mass_list.append(myM1)
                    secondary_mass_list.append(myq * myM1)
                    porb_list.append(10**mylogP)
                    ecc_list.append(mye)
                    binfrac_list.append(mybinfrac)
                mass_binaries += myM1
                mass_binaries += myq * myM1
                n_binaries += 1
            else:
                single_mass_list.append(myM1)
                mass_singles += myM1
                n_singles += 1

        return

def idl_tabulate(x, f, p=5):
    """Function that replicates the IDL int_tabulated function
    which performs a p-point integration on a tabulated set of data

    Parameters
    ----------
    x : array
        tabulated x-value data
    f : array
        tabulated f-value data, same size as x
    p : int
        number of chunks to divide tabulated data into
        Default: 5

    Returns
    -------
    ret : float
        Integration result
    """

    def newton_cotes(x, f):
        if x.shape[0] < 2:
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = scipy.integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)

    ret = 0
    for idx in range(0, x.shape[0], p - 1):
        ret += newton_cotes(x[idx: idx + p], f[idx: idx + p])
    return ret