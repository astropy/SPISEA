##################################################
#
# Original code was taken from libimf package written by Jan Pflamm-Altenburg
# and has been modified only marginally. The libimf code was licensed under
# a GNU General Public License.
# 
# When I use this code, I should cite Pflamm-Altenburg & Kroupa 2006
#
# Unfortunately, the code was almost completely un-commented, so all
# comments are mine. I have also added substantially to the code to
# make more convenient and general purpose functions.
#
##################################################

import numpy as np
import pdb
import logging

log = logging.getLogger('imf')

defaultMFamp = 0.44
defaultMFindex = 0.51
defaultCSFamp = 0.50
defaultCSFindex = 0.45
defaultCSFmax = 3
multiplicity = 0

def sample_imf(massLimits, imfSlopes, totalMass,
               makeMultiples=True, 
               multiMFamp=defaultMFamp, multiMFindex=defaultMFindex,
               multiCSFamp=defaultCSFamp, multiCSFindex=defaultCSFindex,
               multiCSFmax=defaultCSFmax,
               multiQindex=-0.4, multiQmin=0.01,
               verbose=False):
    """
    Randomly sample from an multi-part powerlaw IMF with specified slope mass
    limits until the desired total mass is reached. The maximum
    stellar mass is not allowed to exceed the total cluster mass.
    The simulated total mass will not be exactly equivalent to the
    desired total mass; but we will take one star above or below
    (whichever brings us closer to the desired total) the desired
    total cluster mass break point.

    IMF Slope is -2.35 for Salpeter.

    To specify an multi-part powerlaw with N segments:
    massLimits - numpy array of size N+1 with the upper and lower mass limits for each segment
    imfSlopes - numpy array of size N with the power law slopes (alpha) for each segment.
    """

    if (massLimits[-1] > totalMass) and verbose:
        print 'sample_imf: Setting maximum allowed mass to %d' % \
            (totalMass)

        massLimits[-1] = totalMass

    imf = IMF_broken_powerlaw(massLimits, imfSlopes)
    imf.imf_norm_cl_wk04(totalMass)

    # First estimate the mean number of stars expected
    meanNumber = imf.imf_int_xi(massLimits[0], massLimits[-1])

    simTotalMass = 0
    newStarCount = round(meanNumber)
    if not makeMultiples:
        newStarCount *= 1.1

    masses = np.array([], dtype=float)
    isMultiple = np.array([], dtype=bool)
    compMasses = []
    systemMasses = np.array([], dtype=float)

    loopCnt = 0

    while simTotalMass < totalMass:
        # Generate a random distribution 20% larger than
        # the number we expect to need.
        uniX = np.random.rand(newStarCount)

        # Convert into the IMF from the inverted CDF
        newMasses = imf.imf_dice_star_cl(uniX)

        if makeMultiples:
            compMasses = [[] for ii in range(len(newMasses))]

            # Determine the multiplicity of every star
            MF, CSF = binary_properties(newMasses, MFamp=multiMFamp, MFindex=multiMFindex,
                                        CSFamp=multiCSFamp, CSFindex=multiCSFindex, CSFmax=multiCSFmax)
            newIsMultiple = np.random.rand(newStarCount) < MF
            newSystemMasses = newMasses.copy()
        
            # Calculate number and masses of companions
            for ii in range(len(newMasses)):
                if newIsMultiple[ii]:
                    n_comp = 1 + np.random.poisson((CSF[ii]/MF[ii]) - 1)
                    q_values = q_cdf_inv(np.random.rand(n_comp), multiQmin, multiQindex)
                    m_comp = q_values * newMasses[ii]

                    # Only keep companions that are more than the minimum mass
                    mdx = np.where(m_comp >= massLimits[0])
                    compMasses[ii] = m_comp[mdx]
                    newSystemMasses[ii] += compMasses[ii].sum()

                    # Double check for the case when we drop all companions.
                    # This happens a lot near the minimum allowed mass.
                    if len(mdx) == 0:
                        newIsMultiple[ii] == False

            newSimTotalMass = newSystemMasses.sum()
            isMultiple = np.append(isMultiple, newIsMultiple)
            systemMasses = np.append(systemMasses, newSystemMasses)
        else:
            newSimTotalMass = newMasses.sum()

        # Append to our primary masses array
        masses = np.append(masses, newMasses)

        if (loopCnt >= 0) and verbose:
            print 'sample_imf: Loop %d added %.2e Msun to previous total of %.2e Msun' % \
                (loopCnt, newSimTotalMass, simTotalMass)

        simTotalMass += newSimTotalMass
        newStarCount = meanNumber * 0.1  # increase by 20% each pass
        loopCnt += 1
        
    # Make a running sum of the system masses
    if makeMultiples:
        massCumSum = systemMasses.cumsum()
    else:
        massCumSum = masses.cumsum()

    # Find the index where we are closest to the desired
    # total mass.
    idx = np.abs(massCumSum - totalMass).argmin()

    masses = masses[:idx+1]

    if makeMultiples:
        systemMasses = systemMasses[:idx+1]
        isMultiple = isMultiple[:idx+1]
        compMasses = compMasses[:idx+1]
    else:
        isMultiple = np.zeros(len(masses), dtype=bool)
        systemMasses = masses

    return (masses, isMultiple, compMasses, systemMasses)



class IMF(object):
    def __init__(self, massLimits=np.array([0.1,150]), multiplicity=None):
        """
        The IMF base class. The multiplicity implementation is here.
        """
        self.multi_props = multiplicity
        self.mass_limits = massLimits

    def generateCluster(self, totalMass):
        """
        Generate a cluster of stellar systems with the specified IMF.
        
        Randomly sample from an IMF with specified mass
        limits until the desired total mass is reached. The maximum
        stellar mass is not allowed to exceed the total cluster mass.
        The simulated total mass will not be exactly equivalent to the
        desired total mass; but we will take one star above or below
        (whichever brings us closer to the desired total) the desired
        total cluster mass break point.

        Primary stars are sampled from the IMF, companions are generated
        based on the multiplicity properties provided.

        @param totalMass The total mass of the cluster (including companions).
        """

        if (massLimits[-1] > totalMass):
            log.info('sample_imf: Setting maximum allowed mass to %d' %
                      (totalMass))

        massLimits[-1] = totalMass

        # Estimate the mean number of stars expected.
        self.normalize(totalMass)
        mean_number = self.getProbabilityBetween(massLimits[0], massLimits[-1])
        newStarCount = round(meanNumber)
        if multiplicity == None:
            newStarCount *= 1.1

        # Generate output arrays.
        masses = np.array([], dtype=float)
        isMultiple = np.array([], dtype=bool)
        compMasses = []
        systemMasses = np.array([], dtype=float)

        # Loop through and add stars to the cluster until we get to
        # the desired total cluster mass.
        totalMassTally = 0
        loopCnt = 0

        while totalMassTally < totalMass:
            # Generate a random number array.
            uniX = np.random.rand(newStarCount)

            # Convert into the IMF from the inverted CDF
            newMasses = imf.imf_dice_star_cl(uniX)

            if multiplicity:
                compMasses = [[] for newMass in newMasses]

                # Determine the multiplicity of every star
                MF = multiplicity.getMultiplicityFraction()
                CSF = multiplicity.getCompanionStarFraction()
                
                newIsMultiple = np.random.rand(newStarCount) < MF

                # Copy over the primary masses. Eventually add the companions.
                newSystemMasses = newMasses.copy()
        
                # Calculate number and masses of companions
                for ii in range(len(newMasses)):
                    if newIsMultiple[ii]:
                        # determine number of companions
                        n_comp = 1 + np.random.poisson((CSF[ii]/MF[ii]) - 1)

                        # Determine the mass ratios of the companions
                        q_values = multiplicity.getMassRatios(np.random.rand(n_comp))

                        # Determine the masses of the companions
                        m_comp = q_values * newMasses[ii]

                        # Add in seperation information

                        # Only keep companions that are more than the minimum mass
                        compMasses[ii] = m_comp[m_comp >= massLimits[0]]
                        newSystemMasses[ii] += compMasses[ii].sum()

                        # Double check for the case when we drop all companions.
                        # This happens a lot near the minimum allowed mass.
                        if len(compMasses) == 0:
                            newIsMultiple[ii] == False

                newTotalMassTally = newSystemMasses.sum()
                isMultiple = np.append(isMultiple, newIsMultiple)
                systemMasses = np.append(systemMasses, newSystemMasses)
            else:
                newTotalMassTally = newMasses.sum()

            # Append to our primary masses array
            masses = np.append(masses, newMasses)
            
            if (loopCnt >= 0):
                log.info('sample_imf: Loop %d added %.2e Msun to previous total of %.2e Msun' %
                         (loopCnt, newTotalMassTally, totalMassTally))

            totalMassTally += newTotalMassTally
            newStarCount = meanNumber * 0.1  # increase by 20% each pass
            loopCnt += 1
        
        # Make a running sum of the system masses
        if multiplicity:
            massCumSum = systemMasses.cumsum()
        else:
            massCumSum = masses.cumsum()

        # Find the index where we are closest to the desired
        # total mass.
        idx = np.abs(massCumSum - totalMass).argmin()

        masses = masses[:idx+1]

        if multiplicity:
            systemMasses = systemMasses[:idx+1]
            isMultiple = isMultiple[:idx+1]
            compMasses = compMasses[:idx+1]
        else:
            isMultiple = np.zeros(len(masses), dtype=bool)
            systemMasses = masses

        return (masses, isMultiple, compMasses, systemMasses)
        

        
    
class IMF_broken_powerlaw(IMF):
    def __init__(self, mass_limits, powers, multiplicity=None):
        """
        Initialze a multi-part power-law with N parts. Each part of the
        power-law is described with a probability density function:

            `P(m) \propto m ** power[n]`  

        for mass_limits[n] < m <= mass_limits[n+1].

        Parameters
        ----------
        mass_limits : numpy array
            Array of length (N + 1) with lower and upper limits of 
            the power-law segments.

        coefficients : numpy array
            Array of length N that contains the powers for each
            power-law segment.
        """
        if len(mass_limits) != len(powers) + 1:
            msg = 'Incorrect specification of multi-part powerlaw.\n'
            msg += '    len(massLimts) != len(powers)+1\n'
            msg += '    len(massLimits) = \n' + len(massLimits)
            msg += '    len(powers) = \n' + len(powers)

            raise RuntimeException(msg)

        self._mass_limits = np.atleast_1d(mass_limits)
        self._m_limits_low = mass_limits[0:-1]
        self._m_limits_high = mass_limits[1:]
        self._powers = powers


        # Calculate the coeffs to make the function continuous
        nterms = len(self._powers)
        coeffs = np.ones(nterms, dtype=float)

        # First term is just 1.0
        # Subsequent terms are products of previous terms and then some.
        for i in range(1, nterms):
            y = self.mLimitsLow[i]**powers[i-1]
            z = self.mLimitsLow[i]**powers[i]

            coeffs[i] *= coeffs[i-1] * y / z

        self.nterms = nterms
        self.coeffs = coeffs
        self.k = 1

    def imf_xi(self, m):
        """
        Probability density describing the IMF.

        Input:
        m - mass of a star

        Output:
        xi - probability of measuring that mass.
        """
        returnFloat = type(m) == float
        
        m = np.atleast_1d(m)
        xi = np.zeros(len(m), dtype=float)
        
        for i in range(len(xi)):
            tmp = gamma_closed(m[i], self.mLimitsLow, self.mLimitsHigh)
            tmp *= self.coeffs * m[i]**self.powers
            y = tmp.sum()
            z = delta(m[i] - self.mLimitsHigh).prod()
            xi[i] = self.k * z * y

        if returnFloat:
            return xi[0]
        else:
            return xi

    def imf_mxi(self, m):
        """
        Mass-weighted probability m*xi
        """
        returnFloat = type(m) == float
        m = np.atleast_1d(m)
        mxi = np.zeros(len(m), dtype=float)
        
        for i in range(len(mxi)):
            tmp = gamma_closed(m[i], self.mLimitsLow, self.mLimitsHigh)
            tmp *= self.coeffs * m[i]**(self.powers+1)
            y = tmp.sum()
            z = delta(m[i] - self.mLimitsHigh).prod()
            mxi[i] = self.k * z * y

        if returnFloat:
            return mxi[0]
        else:
            return mxi


    def getProbabilityBetween(self, massLo, massHi):
        """
        Return the integrated probability between some low and high mass value.
        """
        return self.prim_xi(massHi) - self.prim_xi(massLo)
    
    def imf_int_xi(self, left, right):
        """
        Return the integrated probability between some low and high mass value.
        """
        return self.prim_xi(right) - self.prim_xi(left)
    
    def imf_int_mxi(self, left, right):
        """
        Return the total mass between some low and high stellar mass value.
        """
        return self.prim_mxi(right) - self.prim_mxi(left)

    def prim_xi(self, a):
        returnFloat = type(a) == float

        a = np.atleast_1d(a)
        val = np.zeros(len(a), dtype=float)

        for i in range(len(val)):
            t1 = theta_open(a[i] - self._m_limits_high) * self.coeffs
            t2 = imf_prim_power(self._m_limits_high, self._powers)
            t3 = imf_prim_power(self._m_limits_low, self._powers)
            y1 = (t1 * (t2 - t3)).sum()

            t1 = gamma_closed(a[i], self._m_limits_low, self._m_limits_high) 
            t1 *= self.coeffs
            t2 = imf_prim_power(a[i], self._powers)
            t3 = imf_prim_power(self._m_limits_low, self._powers)
            y2 = (t1 * (t2 - t3)).sum()

            val[i] = self.k * (y1 + y2)

        if returnFloat:
            return val[0]
        else:
            return val

    def prim_mxi(self, a):
        returnFloat = type(a) == float
        
        a = np.atleast_1d(a)
        val = np.zeros(len(a), dtype=float)

        for i in range(len(val)):
            t1 = theta_open(a[i] - self._m_limits_high) * self.coeffs
            t2 = imf_prim_power(self._m_limits_high, self._powers+1)
            t3 = imf_prim_power(self._m_limits_low, self._powers+1)
            y1 = (t1 * (t2 - t3)).sum()
            
            t1 = gamma_closed(a[i], self._m_limits_low, self._m_limits_high) 
            t1 *= self.coeffs
            t2 = imf_prim_power(a[i], self._powers+1)
            t3 = imf_prim_power(self._m_limits_low, self._powers+1)
            y2 = (t1 * (t2 - t3)).sum()

            val[i] = self.k * (y1 + y2)

        if returnFloat:
            return val[0]
        else:
            return val

    def normalize(self, Mcl, Mmin=None, Mmax=None):
        """
        Normalize the IMF to a total cluster mass within a specified
        minimum and maximum stellar mass range.
        """
        self.k = 1.0
        self.Mcl = Mcl
        
        if Mmax == None:
            Mmax = self.mLimitsHigh[-1]

        if Mmin == None:
            Mmin = self.mLimitsLow[0]

        if Mmax > Mcl:
            Mmax = Mcl
            
        if Mmax > self.mLimitsHigh[-1]:
            Mmax = self.mLimitsHigh[-1]

        if Mmin < self.mLimitsLow[0]:
            Mmin = self.mLimitsLow[0]

        self.norm_Mmin = Mmin
        self.norm_Mmax = Mmax
        
        self.k = Mcl / self.imf_int_mxi(self.norm_Mmin, self.norm_Mmax)
        self.lamda = self.imf_int_xi_cl(self.mLimitsLow[0], self.massLimits)

    def imf_norm_cl_wk04(self, Mcl, Mmax=None, Mmin=None):
        self.k = 1.0
        self.Mcl = Mcl

        if Mmax == None:
            Mmax = self._m_limits_high[-1]

        if Mmin == None:
            Mmin = self._m_limits_low[0]

        if Mmax > Mcl:
            Mmax = Mcl
            
        if Mmax > self._m_limits_high[-1]:
            Mmax = self._m_limits_high[-1]

        if Mmin < self._m_limits_low[0]:
            Mmin = self._m_limits_low[0]

        a = Mmin
        c = Mmax
        b = (c + a) / 2.0
        while (((c/b)-(a/b)) > 0.00001):
            mb = self.imf_int_mxi(Mmin, b) / self.imf_int_xi(b, Mmax)
            if mb < Mcl:
                a = b
            else:
                c = b
            b = (c + a) / 2.0

        Mmax = b
        self.norm_Mmin = Mmin
        self.norm_Mmax = Mmax

        self.k = Mcl / self.imf_int_mxi(Mmin, Mmax)
        self.lamda = self.imf_int_xi_cl(self._m_limits_low[0], self._mass_limits)

    def imf_xi_cl(self, m):
        return theta_closed(self.norm_Mmax - m) * self.imf_xi(m)

    def imf_mxi_cl(self, m):
        return theta_closed(self.norm_Mmax - m) * self.imf_mxi(m)

    def imf_int_xi_cl(self, left, right):
        t1 = self.prim_xi(right)
        t2 = theta_closed(right - self.norm_Mmax)
        t3 = self.imf_int_xi(self.norm_Mmax, right)
        t4 = self.prim_xi(left)
        t5 = theta_closed(left - self.norm_Mmax)
        t6 = self.imf_int_xi(self.norm_Mmax, left)

        return (t1 - t2*t3) - (t4 - t5*t6)

    def imf_int_mxi_cl(self, left, right):
        t1 = self.prim_mxi(right)
        t2 = theta_closed(right - self.norm_Mmax)
        t3 = self.imf_int_mxi(self.norm_Mmax, right)
        t4 = self.prim_mxi(left)
        t5 = theta_closed(left - self.norm_Mmax)
        t6 = self.imf_int_mxi(self.norm_Mmax, left)

        return (t1 - t2*t3) - (t4 - t5*t6)

    def imf_dice_star_cl(self, r):
        """
        Given a list of random numbers (r), return a list of masses
        selected from the IMF.
        """
        returnFloat = type(r) == float
        r = np.atleast_1d(r)  # Make sure it is an array
        
        x = r * self.lamda[-1]
        y = np.zeros(len(r), dtype=float)
        z = np.ones(len(r), dtype=float)

        for i in range(self.nterms):
            aux = x - self.lamda[i]

            # Only continue for those entries that are in later segments
            idx = np.where(aux >= 0)[0]

            # Maybe we are all done?
            if len(idx) == 0:
                break
            
            x_tmp = x[idx]
            aux_tmp = aux[idx]

            # len(idx) entries
            t1 = aux_tmp / (self.coeffs[i] * self.k)
            t1 += imf_prim_power(self._m_limits_low[i], self._powers[i])
            y_i = gamma_closed(x_tmp, self.lamda[i], self.lamda[i+1])
            y_i *= imf_inv_prim_power(t1, self._powers[i])
            
            # Save results into the y array
            y[idx] += y_i

            z *= delta(x - self.lamda[i])

        if returnFloat:
            return y[0] * z[0]
        else:
            return y * z

class IMFSalpeter1955(IMF_broken_powerlaw):
    def __init__(self, multiplicity=multiplicity):

        massLimits = np.array([0.40, 10.0])
        powers = np.array([-2.3])

        IMF_broken_powerlaw.__init__(self, massLimits, powers,
                                     multiplicity=multiplicity)


class Miller_Scalo_1979(IMF_broken_powerlaw):
    def __init__(self):
        massLimits = np.array([0.1, 1, 10, np.inf])
        powers = np.array([-1.4, -2.5, -3.3])

        IMF_broken_powerlaw.__init__(self, massLimits, powers)

class Kennicutt_1983(IMF_broken_powerlaw):
    def __init__(self):
        massLimits = np.array([0.1, 1, np.inf])
        powers = np.array([-1.4, -2.5])

        IMF_broken_powerlaw.__init__(self, massLimits, powers)

class Kroupa_2001(IMF_broken_powerlaw):
    def __init__(self):
        massLimits = np.array([0.01, 0.08, 0.5, 1, np.inf])
        powers = np.array([-0.3, -1.3, -2.3, -2.3])

        IMF_broken_powerlaw.__init__(self, massLimits, powers)

class Weidner_Kroupa_2004(IMF_broken_powerlaw):
    def __init__(self):
        massLimits = np.array([0.01, 0.08, 0.5, 1, np.inf])
        powers = np.array([-0.3, -1.3, -2.3, -2.35])

        IMF_broken_powerlaw.__init__(self, massLimits, powers)

##################################################
# 
# Generic functions -- see if we can move these up.
#
##################################################
def imf_prim_power(m, power):
    """
    Takes floats or arrays, but returns arrays.
    returns: m**(power + 1) / (power + 1) and handles the case when power = -1
    """
    returnFloat = (type(m) == float) and (type(power) == float)

    m = np.atleast_1d(m)
    power = np.atleast_1d(power)

    z = 1.0 + power
    val = (m**z) / z
    val[power == -1] = np.log(m[power == -1])

    if returnFloat:
        return val[0]
    else:
        return val

def imf_inv_prim_power(x, power):
    """
    returns ((1+power) * x)**(1.0 / (1 + power)) and handles the case 
    when power == -1.
    """
    returnFloat = (type(x) == float) and (type(power) == float)

    x = np.atleast_1d(x)
    power = np.atleast_1d(power)

    z = 1.0 + power
    val = (z * x)**(1.0 / z)
    val[power == -1] = np.exp(x[power == -1])

    if returnFloat:
        return val[0]
    else:
        return val
    

def imf_log_normal(m, mean_logm, sigma_logm):
    returnFloat = (type(m) == float) and (type(mean_logm) == float) and \
        (type(sigma_logm) == float)

    m = np.atleast_1d(m)
    mean_logm = np.atleat_1d(mean_logm)
    sigma_logm = np.atleat_1d(sigma_logm)

    z = np.log10(m) - mean_logm
    val = np.exp(-z**2 / (2.0 * sigma_logm**2)) / m
    
    if returnFloat:
        return val[0]
    else:
        return val

def imf_prim_log_normal(m, mean_logm, sigma_logm):
    returnFloat = (type(m) == float) and (type(mean_logm) == float) and \
        (type(sigma_logm) == float)

    m = np.atleast_1d(m)
    mean_logm = np.atleat_1d(mean_logm)
    sigma_logm = np.atleat_1d(sigma_logm)

    mu = (np.log10(m) - mean_logm) / (1.4142135623731 * sigma_logm)
    val = 2.88586244942136 * sigma_logm * error(mu)
    
    if returnFloat:
        return val[0]
    else:
        return val

def imf_inv_prim_log_normal(x, mean_logm, sigma_logm):
    returnFloat = (type(m) == float) and (type(mean_logm) == float) and \
        (type(sigma_logm) == float)

    m = np.atleast_1d(m)
    mean_logm = np.atleat_1d(mean_logm)
    sigma_logm = np.atleat_1d(sigma_logm)
    
    mu = inv_error(0.346516861952484 * x / sigma_logm)
    val = 10.0**(1.4142135623731 * sigma_logm * mu + mean_logm)
    
    if returnFloat:
        return val[0]
    else:
        return val

def imf_mlog_normal(x, mean_logm, sigma_logm):
    returnFloat = (type(m) == float) and (type(mean_logm) == float) and \
        (type(sigma_logm) == float)

    m = np.atleast_1d(m)
    mean_logm = np.atleat_1d(mean_logm)
    sigma_logm = np.atleat_1d(sigma_logm)

    z = np.log10(m) - mean_logm
    val = np.exp(-z**2 / (2.0 * sigma_logm**2))
    
    if returnFloat:
        return val[0]
    else:
        return val

def imf_prim_mlog_normal(x, mean_logm, sigma_logm):
    returnFloat = (type(m) == float) and (type(mean_logm) == float) and \
        (type(sigma_logm) == float)

    m = np.atleast_1d(m)
    mean_logm = np.atleat_1d(mean_logm)
    sigma_logm = np.atleat_1d(sigma_logm)

    eta = np.log10(m) - mean_logm - (sigma_logm**2 * 2.30258509299405)
    eta /= 1.4142135623731 * sigma_logm

    t1 = (1.15129254649702 * sigma_logm**2) + mean_logm

    val = error(eta)
    val *= 2.88586244942136 * sigma_logm * np.exp(2.30258509299405 * t1)
    
    if returnFloat:
        return val[0]
    else:
        return val
    

def theta_closed(x):
    """
    Pass in float or array of floats (x) and return 0 for x<0
    and 1 for everything else.
    """
    isFloat = type(x) == float

    x = np.atleast_1d(x)
    val = np.ones(len(x), dtype=float)
    val[x < 0] = 0.0

    if isFloat:
        return val[0]
    else:
        return val

def theta_open(x):
    """
    Pass in float or array of floats (x) and return 1 for x>0
    and 0 for everything else.
    """
    isFloat = type(x) == float

    x = np.atleast_1d(x)
    val = np.zeros(len(x), dtype=float)
    val[x > 0] = 1.0

    if isFloat:
        return val[0]
    else:
        return val

def delta(x):
    """
    Pass in float or array of floats (x) and return 0.5 for x==0
    and 1.0 for everything else.
    """
    isFloat = type(x) == float

    x = np.atleast_1d(x)
    val = np.ones(len(x), dtype=float)
    val[x == 0] = 0.5

    if isFloat:
        return val[0]
    else:
        return val

def gamma_closed(m, left, right):
    return theta_closed(m - left) * theta_closed(right - m)


def error(x):
    x2 = x**2
    ax2 = 0.140012288686666 * x2
    
    val = np.sqrt(1.0 - np.exp(-x2*(1.27323954473516+ax2)/(1+ax2)))

    if x >=0:
        return val
    else:
        return -val

def inv_error(x):
    x2 = x**2
    lnx2 = np.log(1.0 - x2)
    aux = 4.54688497944829 + (lnx2 / 2.0)
    y = -aux + np.sqrt(aux**2 - (lnx2 / 0.140012288686666))
    
    val = np.sqrt(y)

    if x>=0:
        return y
    else:
        return -y
    
def binary_properties(mass, MFamp=defaultMFamp, MFindex=defaultMFindex,
                      CSFamp=defaultCSFamp, CSFindex=defaultCSFindex, CSFmax=defaultCSFmax):
    """
    Given a star's mass, determine the probability that the star is in a
    multiple system (multiplicity fraction = MF) and its average number of
    companion stars (companion star fraction = CSF).
    """
    # Multiplicity Fraction
    mf = MFamp * mass**MFindex
    mf[mf > 1] = 1

    # Companion Star Fraction
    csf = CSFamp * mass**CSFindex
    csf[csf > 3] = CSFmax

    return mf, csf

def q_cdf_inv(x, qLo, beta):
    """
    Generative function for companion mass ratio (q = m_comp/m_primary).

    Inputs:
    x -- random number between 0 and 1.
    qLo -- lowest possible mass ratio
    beta -- p(q) \propto q^beta
    """
    b = 1.0 + beta
    return (x * (1.0 - qLo**b) + qLo**b)**(1.0/b)

