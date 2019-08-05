#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Reddening laws.
"""
import pylab as py
import numpy as np
from scipy import interpolate
import pysynphot
from scipy.linalg import solve_banded
import pdb

class RedLawNishiyama09(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Nishiyama et al. 2009 reddening law. The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.5 - 8.0 microns.
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 3-8 microns
        wave = np.arange(0.5, 8.0, 0.001)
        
        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        wave_vals, Alambda_scaled = RedLawNishiyama09.derive_nishiyama09(wave)

        # Convert wavelength to angstrom
        wave_vals *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave_vals, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Nishiyama09',
                                                  litref='Nishiyama+ 2009')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave_vals)
        self.high_lim = max(wave_vals)
        self.name = 'N09'
    
    @staticmethod
    def derive_nishiyama09(wavelength):
        """ 
        Calculate the N09 extinction law as defined in the paper:
        a A_lambda/AKs = power law of exponent -2.0 between JHK. Then
        use a *linear* interpolation in 1/lambda space to go from J to the V-band observation,
        in order to avoid imposing more structure. A cublic spline interpolation
        across the wavelength points is used longward of K-band

        Parameters
        ----------
        wavelength : float
            in microns
        AKs : float
            in magnitudes
        """
        #-----Define power law extinction law between JHK----#
        jhk_idx = np.where( (wavelength >= 1.25) & (wavelength <= 2.14) )
        
        alpha = 2.0
        wave_jhk = wavelength[jhk_idx]

        A_jhk = wave_jhk**(-1.0*alpha)
        A_Ks_jhk = A_jhk / A_jhk[-1]

        #----Now do a linear interpolation (in log(1/lambda) vs log(A/AKs) space) between 1.25 microns and 0.551 microns---#
        jv_idx = np.where( (wavelength < 1.25) & (wavelength > 0.551) )
        Av = 16.13
        func = interpolate.interp1d(np.log10(np.array([1.0/1.25, 1.0/0.551])), np.log10(np.array([A_Ks_jhk[0], Av])),
                                        kind='linear')
        A_Ks_jv = func(np.log10(1.0 / wavelength[jv_idx]))

        # Convert back to linear space
        A_Ks_jv = 10**A_Ks_jv

        #---Do a spline interpolation for the rest of the (long-wavelength) law---#
        # We do this since no other function form is given
        long_idx = np.where(wavelength > 2.14)
        wave = np.array([0.551, 1.25, 1.63, 2.14, 3.545, 4.442, 5.675, 7.760])
        A_AKs = np.array([16.13, 3.02, 1.73, 1.00, 0.500, 0.390, 0.360, 0.430])
        
        spline_interp = interpolate.splrep(wave, A_AKs, k=3, s=0)
        A_AKs_long = interpolate.splev(wavelength[long_idx], spline_interp)
        
        # Stitch together sections for the final law
        wave_vals = np.concatenate((wavelength[jv_idx[0]], wavelength[jhk_idx[0]]))
        A_AKs_vjhk = np.concatenate((A_Ks_jv, A_Ks_jhk))

        # Now add the long-wavelength law
        wave_vals = np.concatenate((wave_vals, wavelength[long_idx[0]]))
        A_AKs_final = np.concatenate((A_AKs_vjhk, A_AKs_long))

        return wave_vals, A_AKs_final

    def Nishiyama09(self, wavelength, AKs):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to derive extinction for, in microns
        AKs : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))
            
        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AKs

        return A_at_wave
        
class RedLawCardelli(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Cardelli et al. 1989 reddening law. The returned object is 
    
    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.3 - 3.0 microns.

    User must specify Rv
    """
    def __init__(self, Rv):
        # Fetch the extinction curve, pre-interpolate across 0.3-3 microns
        wave = np.arange(0.3, 3.0, 0.001)
        
        # This will eventually be scaled by AKs when you
        # call reddening(). Produces A_lambda for AKs = 1, which will be 
        # scaled later. Expects wavelength in microns
        Alambda_scaled = RedLawCardelli.derive_cardelli(wave, Rv)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Cardelli89',
                                                  litref='Cardelli+ 2009')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'C89'

    @staticmethod
    def derive_cardelli(wavelength, Rv):
        """
        Cardelli extinction law. This produces extinction values expected
        for AKs = 1
        """
        x = 1.0 / np.array(wavelength)

        # check for applicability
        if (np.min(x) < 0.3):
            print( 'wavelength is longer than applicable range for Cardelli law')
            return None

        if (np.max(x) > 8.0):
            print( 'wavelength is shorter than applicable range for Cardelli law')
            return None
        
        # Set up some arrays for coefficients that we will need
        a = np.zeros(len(x), dtype=float)
        b = np.zeros(len(x), dtype=float)

        y = x - 1.82

        # Calculate coefficients for long wavelengths (low wavenumber)
        # Wavenumger <= 1.1 (Eq. 2a, 2b)
        idx = np.where(x <= 1.1)[0]
        a[idx] =  0.574 * x[idx] ** 1.61
        b[idx] = -0.527 * x[idx] ** 1.61

        # Calculate coefficients for intermediate wavelengths
        # 1.1 < wavenumber <= 3.3 (Eq. 3a, 3b)
        idx = np.where((x > 1.1) & (x <= 3.3))[0]
        yy = y[idx]
        a[idx] = 1 + (0.17699 * yy) - (0.50447 * yy ** 2) - \
            (0.02427 * yy ** 3) + (0.72085 * yy ** 4) + \
            (0.01979 * yy ** 5) - (0.77530 * yy ** 6) + \
            (0.32999 * yy ** 7)
        b[idx] = (1.41338 * yy) + (2.28305 * yy ** 2) + \
            (1.07233 * yy ** 3) - (5.38434 * yy ** 4) - \
            (0.62251 * yy ** 5) + (5.30260 * yy ** 6) - \
            (2.09002 * yy ** 7)

        # Calculate the long wavelength
        # 3.3 < wavenumber < 5.9 (Eq. 4a, 4b)
        idx = np.where((x > 3.3) & (x < 5.9))[0]
        xx = x[idx]
        a[idx] = 1.752 - (0.316 * xx) - (0.104/((xx - 4.67) ** 2 + 0.341))
        b[idx] = -3.090 + (1.825 * xx) + (1.206/((xx - 4.62) ** 2 + 0.263))

        # Calculate the longest wavelength
        # 5.9 <= wavenumber (Eq. 4a, 4b)
        idx = np.where(x >= 5.9)[0]
        xx = x[idx]
        a[idx] = 1.752 - (0.316 * xx) - (0.104/((xx - 4.67) ** 2 + 0.341)) + \
            (-0.04473 * (xx - 5.9) ** 2) - (0.009779 * (xx - 5.9) ** 3)
        b[idx] = -3.090 + (1.825 * xx) + (1.206/((xx - 4.62) ** 2 + 0.263)) + \
            (0.2130 * (xx - 5.9) ** 2) + (0.1207 * (xx - 5.9) ** 3)

        # A(lam) / A(V), from Eq. 1
        extinction = a + b/Rv

        # Now, want to produce A_lambda / AKs, to match other laws
        k_ind = np.where(abs(x-0.46) == min(abs(x-0.46)))
        Aks_Av = a[k_ind] + b[k_ind]/Rv # Aks / Av
        Av_Aks = 1.0 / Aks_Av # Av / Aks
        
        output = extinction * Av_Aks # (A(lamb) / Av) * (Av / Aks) = (A(lamb) / Aks)

        return output

    def Cardelli89(self, wavelength, AKs):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to derive extinction for, in microns
        AKs : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))
            
        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AKs

        return A_at_wave
    
class RedLawRomanZuniga07(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Roman Zunigoa et al. 2007 reddening law. The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.5 - 8.0 microns.
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.arange(1.0, 8.0, 0.01)
        
        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawRomanZuniga07.derive_romanzuniga07(wave)

        # Convert wavelength to angstrom
        wave *= 10**4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='RomanZuniga07',
                                                  litref='Roman-Zuniga+ 2007')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'RZ07'

    @staticmethod
    def derive_romanzuniga07(wavelength):
        filters = ['J', 'H', 'Ks', '[3.6]', '[4.5]', '[5.8]', '[8.0]']
        wave =      np.array([1.240, 1.664, 2.164, 3.545, 4.442, 5.675, 7.760])
        A_AKs =     np.array([2.299, 1.550, 1.000, 0.618, 0.525, 0.462, 0.455])
        A_AKs_err = np.array([0.530, 0.080, 0.000, 0.077, 0.063, 0.055, 0.059])
        
        # Interpolate over the curve
        spline_interp = interpolate.splrep(wave, A_AKs, k=3, s=0)
        A_AKs_at_wave = interpolate.splev(wavelength, spline_interp)

        return A_AKs_at_wave

    def RomanZuniga07(self, wavelength, AKs):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to derive extinction for, in microns
        AKs : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))
            
        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AKs

        return A_at_wave
    
class RedLawRiekeLebofsky(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Rieke+Lebofsky+85 reddening law. The returned object is 
    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.365 - 13 microns.
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 0.365-13 microns
        wave = np.arange(0.365, 13.0, 0.001)
        
        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawRiekeLebofsky.derive_RiekeLebofsky(wave)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='RiekeLebofsky',
                                                  litref='Rieke+Lebovsky 1985')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'RL85'

    @staticmethod
    def derive_RiekeLebofsky(wavelength):
        """
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Data pulled from Rieke+Lebofsky 1985, Table 3
        """
        filters = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'L', 'M', 
                   '[8.0]', '[8.5]', '[9.0]', '[9.5]', '[10.0]', '[10.5]', 
                   '[11.0]', '[11.5]', '[12.0]', '[12.5]', '[13.0]']
        #wave = np.array([0.365, 0.445, 0.551, 0.658, 0.806, 1.25, 1.635, 2.2, 
        #                 3.77, 4.68, 4.75, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0,
        #                11.5, 12.0, 12.5, 13.0])
        
        # Wavelengths from Nishiyama+09 plot of RL+85 law...slightly different than standard, 
        # drop N filter
        wave = np.array([0.365, 0.445, 0.551, 0.658, 0.806, 1.17, 1.57, 2.12, 
                         3.40, 4.75, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0,
                        11.5, 12.0, 12.5, 13.0])
        A_Av = np.array([1.531, 1.324, 1.00, 0.748, 0.482, 0.282, 0.175, 0.112,
                         0.058, 0.023, 0.02, 0.043, 0.074, 0.087, 0.083,
                         0.074, 0.060, 0.047, 0.037, 0.030, 0.027])
        # Want to change this from A/Av to A/AK
        k_ind = np.where(np.array(filters) == 'K')
        Ak_Av = A_Av[k_ind]
        Av_Ak = 1.0 / Ak_Av

        A_Ak = A_Av * Av_Ak
        
        # Interpolate over the curve
        spline_interp = interpolate.splrep(wave, A_Ak, k=3, s=0)
        A_Ak_at_wave = interpolate.splev(wavelength, spline_interp)

        return A_Ak_at_wave

    def RiekeLebofsky85(self, wavelength, AKs):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to derive extinction for, in microns
        AKs : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))
            
        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AKs

        return A_at_wave

class RedLawDamineli16(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Damineli et al. 2016 reddening law for Wd1. The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.5 - 8.0 microns.
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.arange(0.3, 8.0, 0.001)
        
        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawDamineli16.derive_Damineli16(wave)
        #Alambda_scaled = RedLawDamineli16.derive_Damineli16_old(wave, 1.0)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Damineli16',
                                                  litref='Damineli+ 2016')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'D16'
    

    @staticmethod
    def derive_Damineli16(wavelength):
        """
        Calculate the Damineli+16 extinction law using their equation 19

        Parameters
        ----------
        wavelength : float
            in microns
        AKs : float
            in magnitudes
        """
        # From their eq 19
        x = np.log10(2.159 / wavelength)
        log_A_AKs = -0.015 + 2.33*x + 0.522*x**2. - 3.001*x**3. + 2.034*x**4.

        # Now to convert this back to linear space
        A_AKs_at_wave = 10**log_A_AKs 

        return A_AKs_at_wave

    def Damineli16(self, wavelength, AKs):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to derive extinction for, in microns
        AKs : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))
            
        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AKs

        return A_at_wave

class RedLawDeMarchi16(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    De Marchi et al. 2016 reddening law for 30 Doradus.
    The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.5 - 8.0 microns.
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.arange(0.3, 8.0, 0.001)
        
        # This will eventually be scaled by AK when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawDeMarchi16.derive_DeMarchi16(wave)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='DeMarchi16',
                                                  litref='DeMarchi+ 2016')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'DM16'

    @staticmethod
    def derive_DeMarchi16(wavelength):
        """
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Data pulled from DeMarchi+16, Table 3

        Note: Authors measure R_VI (V) = 3.09 +/- 0.15,
        so we use this to calculate the absolute extinctions in all
        of the other bands. This corresponds to A_I/ A_V = 0.676

        Note that they extrapolate their curve to get to K-band

        Parameters
        ----------
        wavelength : float
            in microns
        AKs : float
            in magnitudes
        """
        AI_AV = 0.676

        # Extracting the values from the paper
        filters = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
        wave = np.array([0.365, 0.445, 0.551, 0.658, 0.806, 1.22, 1.63, 2.19])
        R_VI = np.array([4.41, 3.78, 3.09, 2.58, 2.09, 1.26, 0.84, 0.52])
        R_VI_err = np.array([0.18, 0.15, 0.15, 0.13, 0.17, 0.18, 0.12, 0.08])

        # We'll calculate A_AKs from R_VI
        A_Av = R_VI * (1. - AI_AV)
        AK_Av = A_Av[-1]
        A_AK = A_Av / AK_Av

        # Interpolate over the curve
        spline_interp = interpolate.splrep(wave, A_AK, k=3, s=0)
        A_AK_at_wave = interpolate.splev(wavelength, spline_interp)

        return A_AK_at_wave

    def DeMarchi16(self, wavelength, AK):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to derive extinction for, in microns
        AK : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))
            
        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AK (since law assumes AK = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AK

        return A_at_wave
    
class RedLawFitzpatrick09(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Fitzpactick+09 reddening law.The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.3 - 3.0 microns.

    This reddening law requires 2 free parameters: alpha and R(V).
    Looking at 14 sightlines, the authors generally find alpha ~ 2.5, R(V) ~ 3
    and alpha ~ 1.8, R(V) ~ 5 (Fig 6)
    """
    def __init__(self, alpha, RV):
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.arange(0.7, 3.0, 0.001)
        
        # This will eventually be scaled by AK when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawFitzpactrick09.derive_Fitzpactrick09(wave, alpha, RV)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Fitzpactrick09',
                                                  litref='Fitzpactrick+ 2009')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'F09'

    @staticmethod
    def derive_Fitzpactrick09(wavelength, alpha, RV):
        """
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Data pulled from Fitzpactrick09, equation 5

        Parameters
        ----------
        wavelength : float
            in microns

        alpha: float
            Free parameter alpha

        RV: float
            Free parameter RV
        """
        alpha = float(alpha)
        RV = float(RV)
        
        # First we'll calculate k(lambda - V) = E(lambda - V) / E(B - V),
        # directly from equation 5
        k = (0.349 + 2.087*RV) * (1.0 / (1.0 + (wavelength / 0.507)**alpha)) - RV

        # We'll calculate Alam/Av from K + Rv
        Alam_Av = (k / RV) + 1. 
        
        # Finally, to get A_lambda/Aks we need to divide Alam_Av by AKs_Av.
        # We'll assume central wavelength of 2.14 for Ks
        idx = np.where(abs(wavelength - 2.14) == min(abs(wavelength - 2.14)))

        A_AKs_at_wave = Alam_Av / Alam_Av[idx]

        return A_AKs_at_wave

    def Fitzpactrick09(self, wavelength, AKs):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to derive extinction for, in microns
        AKs : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))
            
        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AKs

        return A_at_wave

class RedLawSchlafly16(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Schafly et al. 2016. The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.5 - 8.0 microns.
    """
    def __init__(self, AH_AKs, x):
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.arange(0.5, 4.8, 0.001)
        
        # This will eventually be scaled by AK when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawSchlafly16.derive_Schlafly16(wave, AH_AKs, x)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Schlafly16',
                                                  litref='Schlafly+ 2016')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'S16'

    @staticmethod
    def derive_Schlafly16(wavelength, AH_AKs, x):
        """
        Calculate Schalfly+16 extinction law according to 
        code provided in appendix of the paper. AH_AKs sets the
        gray component while x sets the shape of the law in an
        Rv-like way
        """
        # Use the function from the Schlafly+16 appendix to get the extinciton law
        # for given AH_AKs and x value. This is given in terms of A_lambda / A(5420)
        law_func = RedLawSchlafly16.Schlafly_appendix(x, AH_AKs)

        # Evaluate function for desired wavelengths (in angstroms)
        law = law_func(wavelength*10**4)
        
        # Now normalize to A_lambda/AKs, rather than A_lambda/A(5420)
        idx = np.where( abs(wavelength - 2.14) == min(abs(wavelength - 2.14)) )
        law_out = law / law[idx]
        
        return law_out

    @staticmethod
    def Schlafly_appendix(x, rhk):
        """ 
        Schlafly+16 extinction law as defined in paper appendix. We've modified
        the wrapper slightly so that the user has control of rhk and x. Here is 
        the comments from that code:
         
        Returns the extinction curve, A(lambda)/A(5420 A), according to
        Schlafly+2016, for the parameter "x," which controls the overall shape of
        the extinction curve in an R(V)-like way.  The extinction curve returned
        is a callable function, which is then invoked with the wavelength, in
        angstroms, of interest.

        The extinction curve is based on broad band photometry between the PS1 g
        band and the WISE W2 band, which have effective wavelengths between 5000
        and 45000 A.  The extinction curve is blindly extrapolated outside that
        range.  The gray component of the extinction curve is fixed by enforcing
        A(H)/A(K) = 1.55 (Indebetouw+2005).  The gray component is relatively
        uncertain, and its variation with x is largely made up.

        Args:
            x: some number controlling the shape of the extinction curve
            ra: extinction vector at anchor wavelengths, default to Schlafly+2016
            dra: derivative of extinction vector at anchor wavelengths, default to
             Schlafly+2016
            lam: anchor wavelengths (angstroms), default to Schlafly+2016

        Returns: the extinction curve E, so the extinction alam = A(lam)/A(5420 A)
        is given by: 
            A = extcurve(x)
            alam = A(lam)
        """
        # Schlafly+2016
        ra = np.array([ 0.65373283,  0.39063843,  0.20197893,  0.07871701, -0.00476316,
                   -0.14213929, -0.23660605, -0.28522577, -0.321301  , -0.33503192])
        dra = np.array([-0.54278669,  0.03404903,  0.36841725,  0.42265873,  0.38247769,
                     0.14148814, -0.04020524, -0.13457319, -0.26883343, -0.36269229])

        # "isoreddening wavelengths" for extinction curve, at E(g-r) = 0.65 reddening
        # T_eff = 4500, Fe/H = 0, log g = 2.5
        lam = np.array([  5032.36441067,   6280.53335141,   7571.85928312,   8690.89321059,
                      9635.52560909,  12377.04268274,  16381.78146718,  21510.20523237,
                     32949.54009328,  44809.4919175 ])


        anchors = ra + x*dra
        # fix gray component so that A(H)/A(K) = 1.55
        anchors += (-anchors[6] + rhk*anchors[7])/(1 - rhk)
        cs0 = CubicSpline(lam, anchors, yp='3d=0')
        # normalize at 5420 angstroms
        return CubicSpline(lam, anchors/cs0(5420.), yp='3d=0')

    def Schlafly16(self, wavelength, AKs):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to derive extinction for, in microns
        AKs : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))
            
        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AKs

        return A_at_wave
        
    @staticmethod
    def derive_schlafly_law(AH):
        """
        Given a value of AH/AK, derive extinction law based on
        Schlafly+16. Color excesses and associated slopes are from
        their Table 3. 

        Problem with this function is that it doesn't take variation in x 
        parameter into account

        Output is an array with Ag/AK, Ar/AK, Ai/AK, Az/AK, Ay/AK, AJ/AK, AH/AK, AKs/AKs=1.0
        """
        # Color excess slopes from Schlafly+16
        slopes = {'gri': 1.395, 'riz': 1.531, 'izy': 1.477, 'zyJ': 0.608 ,
                'yJH': 1.454, 'JHK': 1.943}
        errs = {'gri': 0.014, 'riz': 0.013, 'izy': 0.036, 'zyJ': 0.010,
                'yJH': 0.042, 'JHK': 0.020}

        # Derive the values based on AH
        val = AH - 1
    
        AJ = slopes['JHK'] * val + AH
        Ay = slopes['yJH']*slopes['JHK']*val + AJ
        Az = slopes['zyJ']*slopes['yJH']*slopes['JHK']*val + Ay
        Ai = slopes['izy']*slopes['zyJ']*slopes['yJH']*slopes['JHK']*val + Az
        Ar = slopes['riz']*slopes['izy']*slopes['zyJ']*slopes['yJH']*slopes['JHK']*val + Ai
        Ag = slopes['gri']*slopes['riz']*slopes['izy']*slopes['zyJ']*slopes['yJH']*slopes['JHK']*val + Ar
    
        law = [Ag, Ar, Ai, Az, Ay, AJ, AH, 1.0]

        return law

class RedLawPowerLaw(pysynphot.reddening.CustomRedLaw):
    """
    Power law extinction law, i.e. A_lambda ~ lambda^(-alpha). The
    returned object is a pysynphot CustomRedLaw. NOTE THAT ALPHA INPUT
    IS POSITIVE FOR A NEGATIVE OVERALL EXPONENT

    The wavelength range is 0.5 - 5.0 microns, by default. But user can change
    this with wave_min and wave_max parameters (in microns!)

    Law is normalized such that AKs = 1, but user specifies which
    K-band wavelength to use

    """
    def __init__(self, alpha, K_wave, wave_min=0.5, wave_max=5.0):
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.arange(wave_min, wave_max, 0.001)
        
        # This will eventually be scaled by AK when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawPowerLaw.derive_powerlaw(wave, alpha, K_wave)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Power law')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'pl'

    @staticmethod
    def derive_powerlaw(wavelength, alpha, K_wave):
        """
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Parameters
        ----------
        wavelength : float
            in microns

        alpha: float
            -1.0 * (power law exponent) 
             
        K_wave: float
            Desired K-band wavelength, in microns
        """
        # Create extinction law
        law = wavelength**(-1.0 * alpha)

        # We'll identify K-band as 2.14 microns
        idx = np.where(abs(wavelength - K_wave) == min(abs(wavelength - K_wave)))
        A_AKs_at_wave = law / law[idx]

        return A_AKs_at_wave

    def powerlaw(self, wavelength, AKs):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to derive extinction for, in microns
        AKs : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))
            
        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AKs

        return A_at_wave

class RedLawFritz11(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Fritz+11 reddening law. The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    1.0 - 19 microns.
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 3-8 microns
        wave = np.arange(1.0, 19, 0.001)
        
        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawFritz11.derive_Fritz11(wave)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Nishiyama09',
                                                  litref='Nishiyama+ 2009')
        
        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'F11'
        
    @staticmethod
    def derive_Fritz11(wavelength):
        """
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized to A_Ks = 1

        Data pulled from Fritz+11, Table 2

        Parameters
        ----------
        wavelength : float
            Wavelength range to derive extinction law over, in microns
        """
        # Extinction law definition
        wave = np.array([1.282, 1.736, 2.166, 2.625, 2.758, 2.873, 3.039, 3.297, 3.74, 3.819, 3.907, 4.052,
                             4.376, 5.128, 5.908, 6.772, 7.459, 7.502, 8.76, 12.371, 19.062])
        A_AKs = np.array([7.91, 4.30, 2.49, 1.83, 1.51, 1.84, 2.07, 1.66, 1.19, 1.19, 1.09, 1.01, 1.09, 0.99,
                              1.04, 0.84, 0.81, 0.79, 2.04, 1.34, 1.34])


        # Interpolate over the curve
        spline_interp = interpolate.splrep(wave, A_AKs, k=3, s=0)
        A_at_wave = interpolate.splev(wavelength, spline_interp)

        # We'll call 2.14 microns the K-band
        idx = np.where( abs(wavelength - 2.14) == min(abs(wavelength - 2.14)) )
        A_AKs_at_wave = A_at_wave / A_at_wave[idx] 

        return A_AKs_at_wave

    def Fritz11(self, wavelength, AKs):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelengths to derive extinction for, in microns
        AKs : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))
            
        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AKs

        return A_at_wave
        
class RedLawHosek18(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Hosek+18 reddening law (Wd1 + Arches RC stars). The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.7 - 3.545 microns.
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 3-8 microns
        wave = np.arange(0.7, 3.545, 0.001)
        
        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawHosek18.derive_Hosek18(wave)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Hosek+18',
                                                  litref='Hosek+ 2018')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'H18'
        
    @staticmethod
    def derive_Hosek18(wavelength):
        """ 
        Derive the Hosek+18 extinction law, using the data from Table 4. 
        
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Data pulled from Hosek+18, Table 4

        Parameters
        ----------
        wavelength : float
            Wavelength range to define extinction law over, in microns
        """
        # Extinction law definition
        wave = np.array([0.8059, 0.962, 1.25, 1.53, 2.14, 3.545])
        A_AKs = np.array([9.66, 6.29, 3.56, 2.33, 1.0, 0.50])
        

        # Following Hosek+18, Interpolate over the curve with cubic spline interpolation
        spline_interp = interpolate.splrep(wave, A_AKs, k=3, s=0)
        A_AKs_at_wave = interpolate.splev(wavelength, spline_interp)

        # This curve already assumes A_Ks = 1.0, so we can go straight to
        # output        
        return A_AKs_at_wave

    def Hosek18(self, wavelength, AKs):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to derive extinction for, in microns
        AKs : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))
            
        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AKs

        return A_at_wave

class RedLawHosek18b(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Hosek+18b (Arches IMF) reddening law (Wd1 + Arches RC stars). The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.7 - 3.545 microns.
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 3-8 microns
        wave = np.arange(0.7, 3.545, 0.001)
        
        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawHosek18b.derive_Hosek18b(wave)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Hosek+18b',
                                                  litref='Hosek+ 2018b')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'H18b'
        
    @staticmethod
    def derive_Hosek18b(wavelength):
        """ 
        Derive the Hosek+18 extinction law, using the data from Table 4. 
        
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Data pulled from Hosek+18, Table 4

        Parameters
        ----------
        wavelength : float
            Wavelength range to define extinction law over, in microns
        """
        # Extinction law definition
        wave = np.array([0.8059, 0.962, 1.25, 1.53, 2.14, 3.545])
        A_AKs = np.array([7.943, 5.715, 3.142, 2.04, 1.0, 0.50])
        
        # Following Hosek+18, Interpolate over the curve with cubic spline interpolation
        spline_interp = interpolate.splrep(wave, A_AKs, k=3, s=0)
        A_AKs_at_wave = interpolate.splev(wavelength, spline_interp)

        # This curve already assumes A_Ks = 1.0, so we can go straight to
        # output        
        return A_AKs_at_wave

    def Hosek18b(self, wavelength, AKs):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to derive extinction for, in microns
        AKs : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))
            
        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AKs

        return A_at_wave

class RedLawNoguerasLara18(RedLawPowerLaw):
    """
    An object that returns the reddening vs. wavelength for the 
    Nogueras-Lara+18 reddening law. The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    which consists of a power law with alpha = 2.3

    The wavelength range over which this law is calculated is
    0.8 - 2.8 microns. The law is derived from HAWK-I JHKs 
    observations
    """
    def __init__(self):
        wave_min = 0.8
        wave_max = 2.8
        RedLawPowerLaw.__init__(self, 2.30, 2.15, wave_min=wave_min, wave_max=wave_max)
        
        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = wave_min*10**4
        self.high_lim = wave_max*10**4
        self.name = 'NL18'

    def NoguerasLara18(self, wavelength, AKs):
        """ 
        Return the value of the extinction law at given wavelengths
        with a total overall extinction.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to derive extinction for, in microns
        AKs : float
            in magnitudes
        """
        # If input entry is a single float, turn it into an array
        try:
            len(wavelength)
        except:
            wavelength = [wavelength]

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim*10**-4)) | (max(wavelength) > (self.high_lim*10**-4))):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))    

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * (10**-4)
        law = self.obscuration

        # Find the value of the law at the closest points
        # to wavelength
        A_AKs_at_wave = []
        for ii in wavelength:
            idx = np.where( abs(wave - ii) == min(abs(wave - ii)) )
            A_AKs_at_wave.append(law[idx][0])

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = np.array(A_AKs_at_wave) * AKs

        return A_at_wave    

#---------------------------#
# Cubic spline function from Schalfly+16 appendix
#---------------------------#
def splint(spl, x):
    npts = len(spl.x)
    lo = np.searchsorted(spl.x, x)-1
    lo = np.clip(lo, 0, npts-2)
    hi = lo + 1
    dx = spl.x[hi] - spl.x[lo]
    a = (spl.x[hi] - x)/dx
    b = (x-spl.x[lo])/dx
    y = (a*spl.y[lo]+b*spl.y[hi]+
         ((a**3-a)*spl.y2[lo]+(b**3-b)*spl.y2[hi])*dx**2./6.)
    return y

class CubicSpline:
    def __init__(self, x, y, yp=None):
        npts = len(x)
        mat = np.zeros((3, npts))
        # enforce continuity of 1st derivatives
        mat[1,1:-1] = (x[2:  ]-x[0:-2])/3.
        mat[2,0:-2] = (x[1:-1]-x[0:-2])/6.
        mat[0,2:  ] = (x[2:  ]-x[1:-1])/6.
        bb = np.zeros(npts)
        bb[1:-1] = ((y[2:  ]-y[1:-1])/(x[2:  ]-x[1:-1]) -
                    (y[1:-1]-y[0:-2])/(x[1:-1]-x[0:-2]))
        if yp is None: # natural cubic spline
            mat[1,0] = 1.
            mat[1,-1] = 1.
            bb[0] = 0.
            bb[-1] = 0.
        elif yp == '3d=0':
            mat[1, 0] = -1./(x[1]-x[0])
            mat[0, 1] =  1./(x[1]-x[0])
            mat[1,-1] =  1./(x[-2]-x[-1])
            mat[2,-2] = -1./(x[-2]-x[-1])
            bb[ 0] = 0.
            bb[-1] = 0.
        else:
            mat[1, 0] = -1./3.*(x[1]-x[0])
            mat[0, 1] = -1./6.*(x[1]-x[0])
            mat[2,-2] =  1./6.*(x[-1]-x[-2])
            mat[1,-1] =  1./3.*(x[-1]-x[-2])
            bb[ 0] = yp[0]-1.*(y[ 1]-y[ 0])/(x[ 1]-x[ 0])
            bb[-1] = yp[1]-1.*(y[-1]-y[-2])/(x[-1]-x[-2])
        y2 = solve_banded((1,1), mat, bb)
        self.x, self.y, self.y2 = (x, y, y2)
    def __call__(self, x):
        return splint(self, x)
