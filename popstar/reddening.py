#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Reddening laws.
"""
import pylab as py
import numpy as np
from scipy import interpolate
import pysynphot
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
        wave = np.arange(0.3, 8.0, 0.001)
        
        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawNishiyama09.nishiyama09(wave, 1.0)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Nishiyama09',
                                                  litref='Nishiyama+ 2009')

    @staticmethod
    def nishiyama09(wavelength, AKs, makePlot=False):
        """
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Data pulled from Nishiyama et al. 2009, Table 1

        Parameters
        ----------
        wavelength : float
            in microns
        AKs : float
            in magnitudes
        """
        # Using the SIRIUS filters (MKO), as defined in paper
        filters = ['V', 'J', 'H', 'Ks', '[3.6]', '[4.5]', '[5.8]', '[8.0]']
        wave = np.array([0.551, 1.25, 1.63, 2.14, 3.545, 4.442, 5.675, 7.760])
        A_AKs = np.array([16.13, 3.02, 1.73, 1.00, 0.500, 0.390, 0.360, 0.430])
        A_AKs_err = np.array([0.04,  0.04, 0.03, 0.00, 0.010, 0.010, 0.010, 0.010])

        # Law with HST + VISTA filter (values derived using spline interpolation)
        #filters = ['V', 'F814W', 'Z', 'Y', 'J', 'F160W', 'H', 'Ks', '[3.6]', '[4.5]', '[5.8]', '[8.0]']
        #wave = np.array([0.551, 0.8059, 0.877, 1.02, 1.25, 1.53, 1.645, 2.14, 3.545, 4.442, 5.675, 7.760])
        #A_AKs = np.array([16.13, 8.8707, 7.4337, 5.1866, 3.02, 1.9256, 1.7032, 1.00, 0.500, 0.390, 0.360, 0.430]) 
        
        # Using the 2MASS JHK filters (values derived using spline interpolation)
        #wave = np.array([0.551, 1.24, 1.66, 2.16, 3.545, 4.442, 5.675, 7.760])
        #A_AKs = np.array([16.13, 2.89, 1.62, 1.00, 0.500, 0.390, 0.360, 0.430])        

        # Interpolate over the curve
        spline_interp = interpolate.splrep(wave, A_AKs, k=3, s=0)

        A_AKs_at_wave = interpolate.splev(wavelength, spline_interp)
        A_at_wave = AKs * A_AKs_at_wave

        if makePlot:
            py.clf()
            py.errorbar(wave, A_AKs, yerr=A_AKs_err, fmt='bo', 
                        markerfacecolor='none', markeredgecolor='blue',
                        markeredgewidth=2)

            # Make an interpolated curve.
            wavePlot = np.arange(wave.min(), wave.max(), 0.1)
            extPlot = interpolate.splev(wavePlot, spline_interp)
            py.loglog(wavePlot, extPlot, 'k-')

            # Plot a marker for the computed value.
            py.plot(wavelength, A_AKs_at_wave, 'rs',
                    markerfacecolor='none', markeredgecolor='red',
                    markeredgewidth=2)
            py.xlabel('Wavelength (microns)')
            py.ylabel('Extinction (magnitudes)')
            py.title('Nishiyama et al. 2009')

        return A_at_wave


class RedLawCardelli(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Cardelli et al. 1989 reddening law. The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.5 - 3.0 microns.

    User must specify Rv
    """
    def __init__(self, Rv):
        # Fetch the extinction curve, pre-interpolate across 0.5-3 microns
        wave = np.arange(0.3, 3.0, 0.001)
        
        # This will eventually be scaled by AKs when you
        # call reddening(). Produces A_lambda for AKs = 1, which will be 
        # scaled later. Expects wavelength in microns
        Alambda_scaled = RedLawCardelli.cardelli(wave, Rv, 1.0)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Cardelli89',
                                                  litref='Cardelli+ 2009')

    @staticmethod
    def cardelli(wavelength, Rv, AKs):
        """
        Cardelli extinction law. This produces extinction values expected
        for AKs
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
        # Wavenumger <= 1.1
        idx = np.where(x <= 1.1)[0]
        a[idx] =  0.574 * x[idx] ** 1.61
        b[idx] = -0.527 * x[idx] ** 1.61
        
        # Calculate coefficients for intermediate wavelengths
        # 1.1 < wavenumber <= 3.3
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
        # 3.3 < wavenumber < 5.9
        idx = np.where((x > 3.3) & (x < 5.9))[0]
        xx = x[idx]
        a[idx] = 1.752 - (0.316 * xx) - (0.104/((xx - 4.67) ** 2 + 0.341))
        b[idx] = -3.090 + (1.825 * xx) + (1.206/((xx - 4.62) ** 2 + 0.263))

        # Calculate the longest wavelength
        # 5.9 <= wavenumber
        idx = np.where(x >= 5.9)[0]
        xx = x[idx]
        a[idx] = 1.752 - (0.316 * xx) - (0.104/((xx - 4.67) ** 2 + 0.341)) + \
            (-0.04473 * (xx - 5.9) ** 2) - (0.009779 * (xx - 5.9) ** 3)
        b[idx] = -3.090 + (1.825 * xx) + (1.206/((xx - 4.62) ** 2 + 0.263)) + \
            (0.2130 * (xx - 5.9) ** 2) + (0.1207 * (xx - 5.9) ** 3)

        # A(lam) / A(V)
        extinction = a + b/Rv

        # Now, want to produce A_lambda / AKs, to match other laws
        k_ind = np.where(abs(x-0.46) == min(abs(x-0.46)))
        Aks_Av = a[k_ind] + b[k_ind]/Rv # Aks / Av
        Av_Aks = 1.0 / Aks_Av # Av / Aks
        
        output = extinction * Av_Aks # (A(lamb) / Av) * (Av / Aks) = (A(lamb) / Aks)

        # Up to this point, AKs = 1; now let's scale to user-defined Aks
        output *= AKs

        return output

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
        Alambda_scaled = RedLawRomanZuniga07.romanzuniga07(wave, 1.0, makePlot=False)

        # Convert wavelength to angstrom
        wave *= 10**4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='RomanZuniga07',
                                                  litref='Roman-Zuniga+ 2007')

    @staticmethod
    def romanzuniga07(wavelength, AKs, makePlot=False):
        filters = ['J', 'H', 'Ks', '[3.6]', '[4.5]', '[5.8]', '[8.0]']
        wave =      np.array([1.240, 1.664, 2.164, 3.545, 4.442, 5.675, 7.760])
        A_AKs =     np.array([2.299, 1.550, 1.000, 0.618, 0.525, 0.462, 0.455])
        A_AKs_err = np.array([0.530, 0.080, 0.000, 0.077, 0.063, 0.055, 0.059])
        
        # Interpolate over the curve
        spline_interp = interpolate.splrep(wave, A_AKs, k=3, s=0)

        A_AKs_at_wave = interpolate.splev(wavelength, spline_interp)
        A_at_wave = AKs * A_AKs_at_wave

        if makePlot:
            py.clf()
            py.errorbar(wave, A_AKs, yerr=A_AKs_err, fmt='bo', 
                        markerfacecolor='none', markeredgecolor='blue',
                        markeredgewidth=2)

            # Make an interpolated curve.
            wavePlot = np.arange(wave.min(), wave.max(), 0.1)
            extPlot = interpolate.splev(wavePlot, spline_interp)
            py.loglog(wavePlot, extPlot, 'k-')

            # Plot a marker for the computed value.
            py.plot(wavelength, A_AKs_at_wave, 'rs',
                    markerfacecolor='none', markeredgecolor='red',
                    markeredgewidth=2)
            py.xlabel('Wavelength (microns)')
            py.ylabel('Extinction (magnitudes)')
            py.title('Roman Zuniga et al. 2007')


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
        Alambda_scaled = RedLawRiekeLebofsky.RiekeLebofsky(wave, 1.0)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='RiekeLebofsky',
                                                  litref='Rieke+Lebovsky 1985')

    @staticmethod
    def RiekeLebofsky(wavelength, AKs, makePlot=False):
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
        A_at_wave = AKs * A_Ak_at_wave

        if makePlot:
            py.clf()
            py.plot(wave, A_Ak, 'r.', ms = 8)

            # Make an interpolated curve.
            wavePlot = np.arange(wave.min(), wave.max(), 0.1)
            extPlot = interpolate.splev(wavePlot, spline_interp)
            py.plot(wavePlot, extPlot, 'k-')
            py.xlabel('Wavelength (microns)')
            py.ylabel('Extinction (magnitudes)')
            py.title('Rieke+Lebofsky 1985')

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
        Alambda_scaled = RedLawDamineli16.Damineli16(wave, 1.0)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Damineli16',
                                                  litref='Damineli+ 2016')

    @staticmethod
    def Damineli16(wavelength, AKs, makePlot=False):
        """
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Data pulled from Daminei+16, Table 1 (RC + Wd1)

        Parameters
        ----------
        wavelength : float
            in microns
        AKs : float
            in magnitudes
        """
        # Using the filters as defined in paper
        filters = ['B', 'V', 'R', 'I', 'Z', 'Y', 'J', 'H', 'Ks', 'W1', 'W2']
        wave = np.array([0.442, 0.537, 0.664, 0.805, 0.878, 1.021, 1.244, 1.651, 2.159, 3.295, 4.481])
        A_AKs = np.array([21.43, 14.95, 11.25, 8.72, 7.23, 5.10, 3.23, 1.77, 1.0, 0.39, 0.26])

        # Interpolate over the curve
        spline_interp = interpolate.splrep(wave, A_AKs, k=3, s=0)

        A_AKs_at_wave = interpolate.splev(wavelength, spline_interp)
        A_at_wave = AKs * A_AKs_at_wave

        if makePlot:
            py.clf()
            py.errorbar(wave, A_AKs, yerr=A_AKs_err, fmt='bo', 
                        markerfacecolor='none', markeredgecolor='blue',
                        markeredgewidth=2)

            # Make an interpolated curve.
            wavePlot = np.arange(wave.min(), wave.max(), 0.1)
            extPlot = interpolate.splev(wavePlot, spline_interp)
            py.loglog(wavePlot, extPlot, 'k-')

            # Plot a marker for the computed value.
            py.plot(wavelength, A_AKs_at_wave, 'rs',
                    markerfacecolor='none', markeredgecolor='red',
                    markeredgewidth=2)
            py.xlabel('Wavelength (microns)')
            py.ylabel('Extinction (magnitudes)')
            py.title('Damineli et al. 2016')

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
        Alambda_scaled = RedLawDeMarchi16.DeMarchi16(wave, 1.0)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='DeMarchi16',
                                                  litref='DeMarchi+ 2016')

    @staticmethod
    def DeMarchi16(wavelength, AK, makePlot=False):
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
        A_at_wave = AK * A_AK_at_wave

        if makePlot:
            py.clf()
            py.errorbar(wave, A_AK, yerr=A_AK_err, fmt='bo', 
                        markerfacecolor='none', markeredgecolor='blue',
                        markeredgewidth=2)

            # Make an interpolated curve.
            wavePlot = np.arange(wave.min(), wave.max(), 0.1)
            extPlot = interpolate.splev(wavePlot, spline_interp)
            py.loglog(wavePlot, extPlot, 'k-')

            # Plot a marker for the computed value.
            py.plot(wavelength, A_AK_at_wave, 'rs',
                    markerfacecolor='none', markeredgecolor='red',
                    markeredgewidth=2)
            py.xlabel('Wavelength (microns)')
            py.ylabel('Extinction (magnitudes)')
            py.title('De Marchi et al. 2016')

        return A_at_wave
    
class RedLawFitzpactrick09(pysynphot.reddening.CustomRedLaw):
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
        wave = np.arange(0.3, 3.0, 0.001)
        
        # This will eventually be scaled by AK when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawFitzpactrick09.Fitzpactrick09(wave, alpha, RV, 1.0)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Fitzpactrick09',
                                                  litref='Fitzpactrick+ 2009')

    @staticmethod
    def Fitzpactrick09(wavelength, alpha, RV, AKs):
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
            
        AKs : float
            in magnitudes
        """
        alpha = float(alpha)
        RV = float(RV)
        
        # First we'll calculate k(lambda - V) = E(lambda - V) / E(B - V),
        # directly from equation 5
        k = (0.349 + 2.087*RV) * (1.0 / (1.0 + (wavelength / 0.507)**alpha)) - RV

        # We'll calculate Alam/Av from K + Rv
        Alam_Av = (k / RV) + 1
        
        # Finally, to get A_lambda/Aks we need to divide Alam_Av by AKs_Av.
        # We'll assume central wavelength of 2.14 for Ks
        idx = np.where(abs(wavelength - 2.14) == min(abs(wavelength - 2.14)))

        A_AKs_at_wave = Alam_Av / Alam_Av[idx]
            
        # Now scale to appropriate overall AKs
        A_at_wave = AKs * A_AKs_at_wave

        return A_at_wave

class RedLawSchlafly16(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Schafly et al. 2016. The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.5 - 8.0 microns.
    """
    def __init__(self, AH_AKs):
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.arange(0.3, 3.0, 0.001)
        
        # This will eventually be scaled by AK when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawSchlafly16.Schlafly16(wave, AH_AKs, 1.0)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Schlafly16',
                                                  litref='Schlafly+ 2016')

    @staticmethod
    def Schlafly16(wavelength, AH_AKs, AKs, makePlot=False):
        """
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Data pulled from Schafly+16, Table 3

        Parameters
        ----------
        wavelength : float
            in microns

        AH_AKs: float
            AH / AKs value
            
        AKs : float
            in magnitudes
        """
        # Extract Alam_Aks values from slopes in paper
        A_AKs = RedLawSchlafly16.derive_schlafly_law(AH_AKs)

        # Associated filters and wavelengths
        filters=['g', 'r', 'i', 'z', 'y', 'J', 'H', 'K']
        wave = np.array([0.5032, 0.6281, 0.7572, 0.8691, 0.9636, 1.2377, 1.6382, 2.1510])
        
        # Interpolate over the curve
        spline_interp = interpolate.splrep(wave, A_AKs, k=3, s=0)

        A_AKs_at_wave = interpolate.splev(wavelength, spline_interp)
        A_at_wave = AKs * A_AKs_at_wave

        if makePlot:
            py.clf()
            py.errorbar(wave, A_AKs, yerr=A_AK_err, fmt='bo', 
                        markerfacecolor='none', markeredgecolor='blue',
                        markeredgewidth=2)

            # Make an interpolated curve.
            wavePlot = np.arange(wave.min(), wave.max(), 0.1)
            extPlot = interpolate.splev(wavePlot, spline_interp)
            py.loglog(wavePlot, extPlot, 'k-')

            # Plot a marker for the computed value.
            py.plot(wavelength, A_AKs_at_wave, 'rs',
                    markerfacecolor='none', markeredgecolor='red',
                    markeredgewidth=2)
            py.xlabel('Wavelength (microns)')
            py.ylabel('Extinction (magnitudes)')
            py.title('Schalfly et al. 2016')

        return A_at_wave

    @staticmethod
    def derive_schlafly_law(AH):
        """
        Given a value of AH/AK, derive extinction law based on
        Schlafly+16. Color excesses and associated slopes are from
        their Table 3.

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

class RedLawWesterlund1_byhand(pysynphot.reddening.CustomRedLaw):
    """
    An object that represents the reddening vs. wavelength for the 
    Westerlund 1 Extinction law fitted by hand. The returned object is 

    pysynphot.reddenining.CustomRedLaw (ArraySpectralElement)

    The wavelength range over which this law is calculated is
    0.5 - 8.0 microns.

    NOTE: This was derived using outdated zeropoints. For updated
    zeropoints, use RedLawWesterlund1
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.arange(0.5, 8.0, 0.001)
        
        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawWesterlund1.westerlund1(wave, 1.0)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Westerlund1',
                                                  litref='Lu+ 2015')

    @staticmethod
    def westerlund1(wavelength, AKs, makePlot=False):
        """
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Data pulled from Nishiyama et al. 2009, Table 1
        """

        filters = ['V', 'F814W', 'Z', 'Y', 'J', 'F160W', 'H', 'Ks', '[3.6]', '[4.5]', '[5.8]', '[8.0]']

        # Based on Nishiyama+ 2009, but then by-eye fitting to Wd 1 data for 
        # 0.551, 0.814, 1.25 micron.
        #wave = np.array([0.551, 0.814, 1.25, 1.63, 2.14, 3.545, 4.442, 5.675, 7.760])
        #A_AKs = np.array([16.13, 8.15, 3.19, 1.89, 1.00, 0.500, 0.390, 0.360, 0.430])

        # With HST + VISTA Filters
        wave = np.array([0.551, 0.8059, 0.877, 1.02, 1.25, 1.53, 1.645, 2.14, 3.545, 4.442, 5.675, 7.760])
        A_AKs = np.array([16.13, 8.52, 6.0, 4.32, 3.02, 2.07, 1.82, 1.00, 0.500, 0.390, 0.360, 0.430])        

        # Interpolate over the curve
        spline_interp = interpolate.splrep(wave, A_AKs, k=3, s=0)

        A_AKs_at_wave = interpolate.splev(wavelength, spline_interp)
        A_at_wave = AKs * A_AKs_at_wave

        if makePlot:
            py.clf()
            py.errorbar(wave, A_AKs, yerr=A_AKs_err, fmt='bo', 
                        markerfacecolor='none', markeredgecolor='blue',
                        markeredgewidth=2)

            # Make an interpolated curve.
            wavePlot = np.arange(wave.min(), wave.max(), 0.1)
            extPlot = interpolate.splev(wavePlot, spline_interp)
            py.loglog(wavePlot, extPlot, 'k-')

            # Plot a marker for the computed value.
            py.plot(wavelength, A_AKs_at_wave, 'rs',
                    markerfacecolor='none', markeredgecolor='red',
                    markeredgewidth=2)
            py.xlabel('Wavelength (microns)')
            py.ylabel('Extinction (magnitudes)')
            py.title('Lu et al. 2015')


        return A_at_wave
