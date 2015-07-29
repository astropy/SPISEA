import numpy as np
import pylab as py
from popstar import reddening
from popstar import evolution
from popstar import atmospheres as atm
from scipy import interpolate
from scipy import stats
from gcwork import objects
from pysynphot import spectrum
from pysynphot import ObsBandpass
from pysynphot import observation as obs
import pysynphot
from astropy import constants, units
from astropy.table import Table, Column
# START Need to remove these jlu dependencies #
from jlu.util import plfit
from jlu.nirc2 import synthetic as nirc2syn
# STOP Need to remove these jlu dependencies #
import pickle
import time, datetime
import math
import os, glob
import tempfile
import scipy
import matplotlib
import time
import pdb

defaultAKs = 2.4
defaultDist = 8000

def Vega():
    # Use Vega as our zeropoint... assume V=0.03 mag and all colors = 0.0
    vega = atm.get_kurucz_atmosphere(temperature=9550, 
                                     gravity=3.95,
                                     metallicity=-0.5)

    vega = spectrum.trimSpectrum(vega, 8000, 50000)

    # This is (R/d)**2 as reported by Girardi et al. 2002, page 198, col 1.
    # and is used to convert to flux observed at Earth.
    vega *= 6.247e-17 
    
    return vega

vega = Vega()
redlaw = reddening.RedLawNishiyama09()


def make_observed_isochrone_hst(logAge, AKs=defaultAKs, distance=defaultDist, 
                                iso_dir='./',
                                verbose=False,
                                massSampling=10, 
                                filters={'127m': 'wfc3,ir,f127m',
                                         '139m': 'wfc3,ir,f127m',
                                         '153m': 'wfc3,ir,f153m',
                                         'J': 'nirc2,J',
                                         'H': 'nirc2,H',
                                         'K': 'nirc2,K',
                                         'Kp': 'nirc2,Kp',
                                         'L': 'nirc2,Lp',
                                         '814w': 'wfc3,uvis1,f814w',
                                         '125w': 'wfc3,ir,f125w',
                                         '160w': 'wfc3,ir,f160w'}):
    """
    massSampling - Sample the raw isochrone every ## steps. The default
                   is massSampling = 10, which takes every 10th point.
                   The isochrones are already very finely sampled. Must be
                   an integer value.
    """
    startTime = time.time()

    print 'Making isochrone: log(t) = %.2f  AKs = %.2f  dist = %d' % \
            (logAge, AKs, distance)
    print '     Starting at: ', datetime.datetime.now()
    print '     Usually takes ~5 minutes'

    # Define directory where hst_isochrones are made
    outFileFmt = '{0}iso_{1:.2f}_hst_{2:4.2f}_{3:4s}.fits'
    outFile = outFileFmt.format(iso_dir, logAge, AKs, str(distance).zfill(4))

    c = constants

    # Get solar metallicity models for a population at a specific age.
    evol_model = evolution.MergedPisaEkstromParsec()
    evol = evol_model.isochrone(age=10**logAge)  # solar metallicity
    if verbose:
        print 'Elapsed time while getting merged isochrone: ', \
          time.time() - startTime

    #Eliminate cases where log g is less than 0
    idx = np.where(evol['logg'] > 0)
    evol = evol[idx]
    
    # Trim down the table by selecting every Nth point where
    # N = mass sampling factor.
    evol = evol[::massSampling]

    # Determine which stars are WR stars.
    evol['isWR'] = evol['logT'] != evol['logT_WR']
    
    # Setup output arrays, filter functions, and reddening laws
    # for each filter requested.
    filt_list = {}
    red_list = {}
    
    nrows = len(evol)
    for filt_name, filt_str in filters.iteritems():
        # Setup the final output array.
        mag_col = Column(np.zeros(nrows, dtype=float), name='mag'+filt_name)
        evol.add_column(mag_col)

        # Get the filter transmission function
        filt = get_filter_info(filt_str)
        filt_list[filt_name] = filt

        # Make reddening
        red = redlaw.reddening(AKs).resample(filt.wave)
        red_list[filt_name] = red

    # Convert luminosity to erg/s
    L_all = 10**evol['logL'] * c.L_sun # luminsoity in erg/s

    temp = 10**evol['logT'] * units.K

    # Calculate radius
    R_all = np.sqrt(L_all / (4.0 * math.pi * c.sigma_sb * temp**4))

    # For each temperature extract the synthetic photometry.
    for ii in range(len(temp)):
        t2 = time.time()
        gravity = float( evol['logg'][ii] )
        L = float( L_all[ii] / (units.erg / units.s)) # in erg/s
        T = float(  temp[ii] / units.K)               # in Kelvin
        R = float( R_all[ii] / units.pc)              # in pc

        # Get the atmosphere model now. Wavelength is in Angstroms
        star = atm.get_merged_atmosphere(temperature=T, 
                                         gravity=gravity)
        # Trim wavelength range down to JHKL range (0.5 - 4.25 microns)
        star = spectrum.trimSpectrum(star, 5000, 42500)

        # Convert into flux observed at Earth (unreddened)
        star *= (R / distance)**2  # in erg s^-1 cm^-2 A^-1

        # ----------
        # Now to the filter integrations
        # ----------
        t0 = time.time()
        for filt_name, filt_str in filters.iteritems():
            col_name = 'mag' + filt_name
            filt = filt_list[filt_name]
            red = red_list[filt_name]
            evol[col_name][ii] = mag_in_filter(star, filt, red)
        t1 = time.time()
            
        if verbose:
            print 'It took {0:f} for reddening calc'.format(t1 - t0)
            
            fmt = 'M = %7.3f Msun  T = %5d K  R = %2.1f Rsun  logg = %4.2f  '
            vals = [evol['mass'][ii], T, float(R_all[ii] / units.R_sun), evol['logg'][ii]]

            for filt_name, filt_str in filters.iteritems():
                fmt += filt_name + ' = %4.2f  '
                vals.append(evol['mag' + filt_name][ii])
                
            fmt += 'elapsed time = %4s'
            vals.append(time.time() - startTime)
            
            print  fmt % tuple(vals)
                
                
            t3 = time.time()
            print 'It took {0:.2f} s to process one star'.format(t3 - t2)


    evol.write(outFile, overwrite=True)

    if verbose:
        endTime = time.time()
        print '      Time taken: %d seconds' % (endTime - startTime)

    return

def load_isochrone(logAge=6.78, AKs=defaultAKs, distance=defaultDist,
                   iso_dir='./'):
    """
    Wrapper code that loads an hst isochrone or make a new one if it doesn't
    already exist.
    """
    # Define directory where hst_isochrones exist
    inFileFmt = '{0}iso_{1:.2f}_hst_{2:4.2f}_{3:4s}.fits'
    inFile = inFileFmt.format(iso_dir, logAge, AKs, str(distance).zfill(4))

    if not os.path.exists(inFile):
        make_observed_isochrone_hst(logAge=logAge, AKs=AKs, distance=distance,
                                    iso_dir=iso_dir)

    iso = Table.read(inFile)

    return iso

def make_isochrone_grid():
    """
    Helper routine to make a isochrone grid. logAge is
    hardcoded.
    """
    logAge_arr = np.arange(6.0, 6.7, 0.01)

    for i in logAge_arr:
        load_isochrone(i)

    return


# Little helper utility to get all the bandpass/zeropoint info.
def get_filter_info(name, vega=vega):
    if name.startswith('nirc2'):
        tmp = name.split(',')
        filterName = tmp[-1]
        filter = nirc2syn.filters[filterName]
        flux0 = nirc2syn.filter_flux0[filterName]
        mag0 = nirc2syn.filter_mag0[filterName]
    else:
        filter = ObsBandpass(name)

    vega_obs = obs.Observation(vega, filter, binset=filter.wave, force='taper')
    vega_flux = vega_obs.binflux.sum()
    vega_mag = 0.03

    filter.flux0 = vega_flux
    filter.mag0 = vega_mag
    
    return filter

# Little helper utility to get the magnitude of an object through a filter.
def mag_in_filter(star, filter, extinction):
    """
    Assumes that extinction is already resampled to same wavelengths
    as filter.
    """
    star_in_filter = obs.Observation(star, filter*extinction,
                                     binset=filter.wave, force='taper')
    star_flux = star_in_filter.binflux.sum()
    star_mag = -2.5 * math.log10(star_flux / filter.flux0) + filter.mag0

    return star_mag


