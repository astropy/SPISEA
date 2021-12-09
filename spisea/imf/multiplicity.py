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
        else:
            csf[csf > self.CSF_max] = self.CSF_max

        return csf

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
