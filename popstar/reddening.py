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
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.arange(0.5, 8.0, 0.001)
        
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
        """
        filters = ['V', 'J', 'H', 'Ks', '[3.6]', '[4.5]', '[5.8]', '[8.0]']
        wave = np.array([0.551, 1.25, 1.63, 2.14, 3.545, 4.442, 5.675, 7.760])
        A_AKs = np.array([16.13, 3.02, 1.73, 1.00, 0.500, 0.390, 0.360, 0.430])
        A_AKs_err = np.array([0.04,  0.04, 0.03, 0.00, 0.010, 0.010, 0.010, 0.010])

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
    0.5 - 8.0 microns.
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.arange(0.5, 8.0, 0.001)
        
        # This will eventually be scaled by AKs when you
        # call reddening(). Produces A_lambda for AKs = 1, which will be 
        # scaled later. Adopt Rv=3.1
        Alambda_scaled = RedLawCardelli.cardelli(wave, 3.1)

        # Convert wavelength to angstrom
        wave *= 10 ** 4

        pysynphot.reddening.CustomRedLaw.__init__(self, wave=wave, 
                                                  waveunits='angstrom',
                                                  Avscaled=Alambda_scaled,
                                                  name='Nishiyama09',
                                                  litref='Nishiyama+ 2009')
    @staticmethod
    def cardelli(wavelength, Rv):
        """
        Cardelli extinction law. Note: this produces extinction values expected
        for AKs = 1 mag
        """
        x = 1.0 / np.array(wavelength)

        # check for applicability
        if (wavelength.min() < 0.3):
            print 'wavelength is longer than applicable range for Cardelli law'
            return None

        if (wavelength.max() > 8.0):
            print 'wavelength is shorter than applicable range for Cardelli law'
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

        # Now, want to produce redvals at AKs = 1
        k_ind = np.where(abs(x-0.47) == min(abs(x-0.47)))
        Aks_Av = a[k_ind] + b[k_ind]/Rv
        # Av / Aks
        Av_Aks = 1.0 / Aks_Av 
        
        output = extinction * Av_Aks # If Aks = 1, Av/Aks = Av
        
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
        Alambda_scaled = extinction.romanzuniga07(wave, 1.0, makePlot=False)

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
    
