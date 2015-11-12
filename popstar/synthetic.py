import numpy as np
import pylab as py
from popstar import reddening
from popstar import evolution
from popstar import atmospheres as atm
from scipy import interpolate
from scipy import stats
#from gcwork import objects
from pysynphot import spectrum
from pysynphot import ObsBandpass
from pysynphot import observation as obs
import pysynphot
from astropy import constants, units
from astropy.table import Table, Column
# START Need to remove these jlu dependencies #
from popstar.imf import imf
#from jlu.imf import imf
from popstar.utils import objects
# from jlu.gc.gcwork import objects
from nirc2 import synthetic as nirc2syn
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
defaultEvModel = 'mergedPisaEkstromParsec'

#-----------------------------------------------------#
# Default settings for model_young_cluster code
defaultMFamp = 0.44
defaultMFindex = 0.51
defaultCSFamp = 0.50
defaultCSFindex = 0.45
defaultCSFmax = 3
#-----------------------------------------------------#

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

def match_model_mass(isoMasses,theMass):
    dm = np.abs(isoMasses - theMass)
    mdx = dm.argmin()

    # Model mass has to be within 2% of the desired mass
    if (dm[mdx] / theMass) > 0.1:
        return None
    else:
        return mdx

def model_young_cluster_object(resolved=False):
    multi = multiplicity.MultiplicityUnresolved()
    imf = imf.Kroupa_2001(multiplicty=multi)
    evo = evolution.MergedPisaEkstromParsec()
    atm_func = atm.get_merged_atmosphere

    log_age = 6.5
    AKs = 1.0
    distance = 8000.0

    if resolved:
        cluster = ResolvedCluster(log_age, AKs, distance, imf, evo, atm_func)
    else:
        cluster = UnresolvedCluster(log_age, AKs, distance, imf, evo, atm_func)

    # Plot the spectrum of the most massive star
    idx = cluster.mass.argmax()
    plt.clf()
    plt.plot(cluster.stars[idx].wave, cluster.stars[idx].flux, 'k.')

    # Plot an integrated spectrum of the whole cluster.
    wave, flux = cluster.get_integrated_spectrum()
    plt.clf()
    plt.plot(wave, flux, 'k.')

    return

class Cluster(object):
    def __init__(self, logAge, AKs, distance,
                 imf, evolution_model, atmosphere_func,
                 verbose=False): 
        """
        Code to model a cluster with user-specified logAge, AKs, and distance.
        Must also specify directory containing the isochrone (made using popstar
        synthetic code).

        Can also specify IMF slope, mass limits, cluster mass, and parameters for
        multiple stars
        """

        c = constants

        logAgeString = '0%d' % (int(logAge * 100))

        # Sample a power-law IMF randomly
        results = imf.sample_imf(massLimits, imfSlopes, clusterMass,
                                 makeMultiples=makeMultiples,
                                 multiMFamp=MFamp, multiMFindex=MFindex,
                                 multiCSFamp=CSFamp, multiCSFindex=CSFindex,
                                 multiCSFmax=CSFmax,
                                 multiQmin=qMin, multiQindex=qIndex)
        mass = results[0] # masses of the stars in the simulated cluster
        isMultiple = results[1]
        compMasses = results[2]
        systemMasses = results[3]

        

    
class ResolvedCluster(Cluster):
    def __init__(self, logAge, AKs, distance,
                 imf, evolution_model, atmosphere_func,
                 verbose=False):

        # this object is going to be handled by Jessica/Matt
        # don't bother fixing anything further in here
        
        Iso = load_isochrone(logAge=logAge, AKs=AKs, distance=distance, evModel=evModel,
                                       iso_dir=iso_dir, mag_calc=True)

        mag814w = np.zeros(len(Cluster.mass), dtype=float)
        mag127m = np.zeros(len(Cluster.mass), dtype=float)
        mag139m = np.zeros(len(Cluster.mass), dtype=float)
        mag153m = np.zeros(len(Cluster.mass), dtype=float)
        magJ = np.zeros(len(Cluster.mass), dtype=float)
        magH = np.zeros(len(Cluster.mass), dtype=float)
        magK = np.zeros(len(Cluster.mass), dtype=float)
        magKp = np.zeros(len(Cluster.mass), dtype=float)
        magL = np.zeros(len(Cluster.mass), dtype=float)
        temp = np.zeros(len(Cluster.mass), dtype=float)
        logg = np.zeros(len(Cluster.mass), dtype=float)
        logL = np.zeros(len(Cluster.mass), dtype=float)
        isWR = np.zeros(len(Cluster.mass), dtype=bool)   

        for ii in range(len(Cluster.mass)):
            # Find the closest model mass (returns None, if nothing with dm = 0.1
            mdx = match_model_mass(Iso.mass,Cluster.mass[ii])
            if mdx == None:
                continue

            # if end up making the magnitudes attributes of Iso in make/load
            # isochrone, update this section as in the first one
            mag814w[ii] = Iso.mag814w[mdx]
            mag127m[ii] = iso['mag127m'][mdx]
            mag139m[ii] = iso['mag139m'][mdx]
            mag153m[ii] = iso['mag153m'][mdx]
            magJ[ii] = iso['magJ'][mdx]
            magH[ii] = iso['magH'][mdx]
            magK[ii] = iso['magK'][mdx]
            magKp[ii] = iso['magKp'][mdx]
            magL[ii] = iso['magL'][mdx]

            temp[ii] = iso['logT'][mdx]
            logg[ii] = iso['logg'][mdx]
            logL[ii] = iso['logL'][mdx]
            isWR[ii] = iso['isWR'][mdx]


            # Determine if this system is a binary.
            if isMultiple[ii]:
                n_stars = len(compMasses[ii])
                for cc in range(n_stars):
                    mdx_cc = match_model_mass(compMasses[ii][cc])
                    if mdx_cc != None:
                        f1 = 10**(-mag[ii]/2.5)
                        f2 = 10**(-iso.mag[mdx_cc]/2.5)
                        mag[ii] = -2.5 * np.log10(f1 + f2)
                    else:
                        print 'Rejected a companion %.2f' % compMasses[ii][cc]


        # Get rid of the bad ones
        idx = np.where(temp != 0)[0]
        cdx = np.where(temp == 0)[0]

        if len(cdx) > 0 and verbose:
            print 'Found %d stars out of mass range: Minimum bad mass = %.1f' % \
                (len(cdx), mass[cdx].min())

        mass = mass[idx]
        mag814w = mag814w[idx]
        mag127m = mag127m[idx]
        mag139m = mag139m[idx]
        mag153m = mag153m[idx]
        magJ = magJ[idx]
        magH = magH[idx]
        magK = magK[idx]
        magKp = magKp[idx]
        magL = magL[idx]


        temp = temp[idx]
        logg = logg[idx]
        logL = logL[idx]
        isWR = isWR[idx]

        isMultiple = isMultiple[idx]
        systemMasses = systemMasses[idx]
        if makeMultiples:
            compMasses = [compMasses[ii] for ii in idx]

        idx_noWR = np.where(isWR == False)[0]

        mag127m_noWR = mag127m[idx_noWR]
        num_WR = len(mag127m) - len(idx_noWR)



    return cluster, iso

            
class UnresolvedCluster(Cluster):
    def __init__(self, logAge, AKs, distance,
                 imf, evolution_model, atmosphere_func,
                 verbose=False): 

        Iso = Isochrone(logAge, AKs, distance, evModel)

        self.spec_list = cluster_spec(Iso)
        self.spec_tot_full = add_all_spectra(self.spec_list)
        self.spec = apply_filter(self.spec_tot_full)
        self.mass_all = Iso.mass_all
        self.log_age = Iso.logAge

    def cluster_spec(isochrone)
        # get the masses from the sampled IMF

        temp = np.zeros(len(Cluster.mass), dtype=float)
        # initialize spec object
        #spec = 
    
        for ii in range(len(Cluster.mass)):
            # Find the closest model mass (returns None, if nothing with dm = 0.1
            mdx = match_model_mass(isochrone.mass,Cluster.mass[ii])
            if mdx == None:
                continue

            temp[ii] = isochrone.T_all[mdx]
            spec[ii] = isochrone.spec_list[mdx]

        # Get rid of the bad ones
        idx = np.where(temp != 0)[0]
        cdx = np.where(temp == 0)[0]

        spec = spec[idx]
        
        return spec

    def add_all_spectra(spec_list)
        # sum all spectra in a list

        

    
    def apply_filter(spectrum, filter)
        # trim/resample spectrum to match a given filter
        
        
 

def evModel_select(evModel):
    evModel_switch = {
        'mergedPisaEkstromParsec': evolution.MergedPisaEkstromParsec(),
        'MergedPisaEkstromParsec_norot': evolution.MergedPisaEkstromParsec_norot(),
        # none of the following models are fully set up in evolution.py yet!!!!
        'Geneva': evolution.GenevaStellarEvolution(),
        'Ekstrom': evolution.EkstromStellarEvolution(),
        'Ekstrom_norot': evolution.EkstromStellarEvolution_norot(),
        'Parsec': evolution.ParsecStellarEvolution(),
        'Pisa': evolution.PisaStellarEvolution()
        }
    pickevModel = evModel_switch.get(evModel)
    return pickevModel
        

class Isochrone(object):
    def __init__(logAge, AKs=defaultAKs, distance=defaultDist,
                                evModel=defaultEvModel,
                                massSampling=2):
    """
    massSampling - Sample the raw isochrone every ## steps. The default
                   is massSampling = 10, which takes every 10th point.
                   The isochrones are already very finely sampled. Must be
                   an integer value.
    """

    c = constants
    
    # Get solar metallicity models for a population at a specific age.
    evol_model = evModel_select(evModel)
    evol = evol_model.isochrone(age=10**logAge)  # solar metallicity

    #Eliminate cases where log g is less than 0
    idx = np.where(evol['logg'] > 0)
    evol = evol[idx]
    
    # Trim down the table by selecting every Nth point where
    # N = mass sampling factor.
    evol = evol[::massSampling]

    # Determine which stars are WR stars.
    evol['isWR'] = evol['logT'] != evol['logT_WR']
  
    # Convert luminosity to erg/s
    self.L_all = 10**evol['logL'] * c.L_sun # luminsoity in erg/s

    self.T_all = 10**evol['logT'] * units.K

    self.mass_all = evol['mass'] # masses in solar masses

    # Calculate radius
    self.R_all = np.sqrt(L_all / (4.0 * math.pi * c.sigma_sb * temp**4))

    # Initialize output for stellar spectra
    self.spec_list = []

    # For each temperature extract the synthetic photometry.
    for ii in range(len(self.T_all)):
        gravity = float( evol['logg'][ii] )
        L = float( self.L_all[ii].cgs / (units.erg / units.s)) # in erg/s
        T = float( self.temp[ii] / units.K)               # in Kelvin
        R = float( self.R_all[ii].to('pc') / units.pc)              # in pc

        # Get the atmosphere model now. Wavelength is in Angstroms
        star = atm.get_merged_atmosphere(temperature=T, 
                                         gravity=gravity)
        # Trim wavelength range down to JHKL range (0.5 - 4.25 microns)
        star = spectrum.trimSpectrum(star, 5000, 42500)

        # Convert into flux observed at Earth (unreddened)
        star *= (R / distance)**2  # in erg s^-1 cm^-2 A^-1
        red = redlaw.reddening(AKs).resample(star.wave)  ## check this
        star *= red
        self.spec_list.append(star)

    evol.meta['RedLaw'] = redlaw.name

    return

def load_isochrone(logAge=6.78, AKs=defaultAKs, distance=defaultDist,
                   evModel=defaultEvModel,
                   iso_dir='./', massSampling=3,
                   mag_calc=False, 
                   filters={'127m': 'wfc3,ir,f127m',
                            '139m': 'wfc3,ir,f127m',
                            '153m': 'wfc3,ir,f153m',
                            'J': 'nirc2,J',
                            'H': 'nirc2,H',
                            'K': 'nirc2,K',
                            'Kp': 'nirc2,Kp',
                            'L': 'nirc2,Lp',
                            '814w': 'acs,wfc1,f814w',
                            '125w': 'wfc3,ir,f125w',
                            '160w': 'wfc3,ir,f160w'}):
    """
    Wrapper code that loads an hst isochrone or make a new one if it doesn't
    already exist.
    """

    # Define directory where hst_isochrones exist
    inFileFmt = '{0}iso_{1:.2f}_hst_{2:4.2f}_{3:4s}.fits'
    inFile = inFileFmt.format(iso_dir, logAge, AKs, str(distance).zfill(4))

    if not os.path.exists(inFile):
        make_observed_isochrone_hst(logAge=logAge, AKs=AKs, distance=distance,
                                    iso_dir=iso_dir, massSampling=massSampling,
                                    filters=filters)

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
        filter = nirc2syn.get_filter_info(filterName, vega=vega)[0]
        # tmp = name.split(',')
        # filterName = tmp[-1]
        # filter = nirc2syn.filters[filterName]
        # flux0 = nirc2syn.filter_flux0[filterName]
        # mag0 = nirc2syn.filter_mag0[filterName]
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

def model_young_cluster_new(logAge, AKs, distance, iso_dir,
                            evModel='mergedPisaEkstromParsec',
                            imfSlopes=np.array([-2.35]), massLimits=np.array([1, 150]),
                            clusterMass=1e4, makeMultiples=False,
                            MFamp=defaultMFamp, MFindex=defaultMFindex,
                            CSFamp=defaultCSFamp, CSFindex=defaultCSFindex,
                            CSFmax=defaultCSFmax,
                            qMin=0.01, qIndex=-0.4, verbose=False):
    """
    Code to model a cluster with user-specified logAge, AKs, and distance.
    Must also specify directory containing the isochrone (made using popstar
    synthetic code).

    Can also specify IMF slope, mass limits, cluster mass, and parameters for
    multiple stars
    """
    c = constants

    logAgeString = '0%d' % (int(logAge * 100))

    iso = load_isochrone(logAge=logAge, AKs=AKs, distance=distance, evModel=evModel,
                                   iso_dir=iso_dir)

    # Sample a power-law IMF randomly
    results = imf.sample_imf(massLimits, imfSlopes, clusterMass,
                             makeMultiples=makeMultiples,
                             multiMFamp=MFamp, multiMFindex=MFindex,
                             multiCSFamp=CSFamp, multiCSFindex=CSFindex,
                             multiCSFmax=CSFmax,
                             multiQmin=qMin, multiQindex=qIndex)
    mass = results[0]
    isMultiple = results[1]
    compMasses = results[2]
    systemMasses = results[3]

    mag814w = np.zeros(len(mass), dtype=float)
    mag127m = np.zeros(len(mass), dtype=float)
    mag139m = np.zeros(len(mass), dtype=float)
    mag153m = np.zeros(len(mass), dtype=float)
    magJ = np.zeros(len(mass), dtype=float)
    magH = np.zeros(len(mass), dtype=float)
    magK = np.zeros(len(mass), dtype=float)
    magKp = np.zeros(len(mass), dtype=float)
    magL = np.zeros(len(mass), dtype=float)
    temp = np.zeros(len(mass), dtype=float)
    logg = np.zeros(len(mass), dtype=float)
    logL = np.zeros(len(mass), dtype=float)
    isWR = np.zeros(len(mass), dtype=bool)
    
    def match_model_mass(theMass):
        dm = np.abs(iso['mass'] - theMass)
        mdx = dm.argmin()

        # Model mass has to be within 2% of the desired mass
        if (dm[mdx] / theMass) > 0.1:
            return None
        else:
            return mdx
        
    for ii in range(len(mass)):
        # Find the closest model mass (returns None, if nothing with dm = 0.1
        mdx = match_model_mass(mass[ii])
        if mdx == None:
            continue

        mag814w[ii] = iso['mag814w'][mdx]
        mag127m[ii] = iso['mag127m'][mdx]
        mag139m[ii] = iso['mag139m'][mdx]
        mag153m[ii] = iso['mag153m'][mdx]
        magJ[ii] = iso['magJ'][mdx]
        magH[ii] = iso['magH'][mdx]
        magK[ii] = iso['magK'][mdx]
        magKp[ii] = iso['magKp'][mdx]
        magL[ii] = iso['magL'][mdx]
        
        temp[ii] = iso['logT'][mdx]
        logg[ii] = iso['logg'][mdx]
        logL[ii] = iso['logL'][mdx]
        isWR[ii] = iso['isWR'][mdx]


        # Determine if this system is a binary.
        if isMultiple[ii]:
            n_stars = len(compMasses[ii])
            for cc in range(n_stars):
                mdx_cc = match_model_mass(compMasses[ii][cc])
                if mdx_cc != None:
                    f1 = 10**(-mag[ii]/2.5)
                    f2 = 10**(-iso.mag[mdx_cc]/2.5)
                    mag[ii] = -2.5 * np.log10(f1 + f2)
                else:
                    print 'Rejected a companion %.2f' % compMasses[ii][cc]
        

    # Get rid of the bad ones
    idx = np.where(temp != 0)[0]
    cdx = np.where(temp == 0)[0]

    if len(cdx) > 0 and verbose:
        print 'Found %d stars out of mass range: Minimum bad mass = %.1f' % \
            (len(cdx), mass[cdx].min())

    mass = mass[idx]
    mag814w = mag814w[idx]
    mag127m = mag127m[idx]
    mag139m = mag139m[idx]
    mag153m = mag153m[idx]
    magJ = magJ[idx]
    magH = magH[idx]
    magK = magK[idx]
    magKp = magKp[idx]
    magL = magL[idx]
    
    
    temp = temp[idx]
    logg = logg[idx]
    logL = logL[idx]
    isWR = isWR[idx]
    
    isMultiple = isMultiple[idx]
    systemMasses = systemMasses[idx]
    if makeMultiples:
        compMasses = [compMasses[ii] for ii in idx]

    idx_noWR = np.where(isWR == False)[0]

    mag127m_noWR = mag127m[idx_noWR]
    num_WR = len(mag127m) - len(idx_noWR)

    cluster = objects.DataHolder()
    cluster.mass = mass
    cluster.Teff = temp
    cluster.logg = logg
    cluster.logL = logL
    cluster.isWR = isWR
    cluster.mag814w = mag814w
    cluster.mag127m = mag127m
    cluster.mag139m = mag139m
    cluster.mag153m = mag153m
    cluster.magJ = magJ
    cluster.magH = magH
    cluster.magK = magK
    cluster.magKp = magKp
    cluster.magL = magL
    cluster.isMultiple = isMultiple
    cluster.compMasses = compMasses
    cluster.systemMasses = systemMasses

    cluster.idx_noWR = idx_noWR
    cluster.mag127m_noWR = mag127m_noWR
    cluster.num_WR = num_WR

    # Summary parameters
    cluster.logAge = logAge
    cluster.AKs = AKs
    cluster.distance = distance
    cluster.massLimits = massLimits
    cluster.imfSlopes = imfSlopes
    cluster.sumIMFmass = clusterMass
    cluster.makeMultiples = makeMultiples
    cluster.MFamp = MFamp
    cluster.MFindex = MFindex
    cluster.CSFamp = CSFamp
    cluster.CSFindex = CSFindex
    cluster.CSFmax = CSFmax
    cluster.qMin = qMin
    cluster.qIndex = qIndex
            
    return cluster, iso
