import numpy as np

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
    `Lu+2013 <https://ui.adsabs.harvard.edu/abs/2013ApJ...764..155L/abstract>`_.

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
    
    """
    def __init__(self, 
                 MF_amp=0.44, MF_power=0.51,
                 CSF_amp=0.50, CSF_power=0.45, CSF_max=3,
                 q_power=-0.4, q_min=0.01):
         
        self.MF_amp = MF_amp
        self.MF_pow = MF_power
        self.CSF_amp = CSF_amp
        self.CSF_pow = CSF_power
        self.CSF_max = CSF_max
        self.q_pow = q_power
        self.q_min = q_min

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

        return n_comp
