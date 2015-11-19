import numpy as np
import pylab as py
from popstar import reddening
from popstar import evolution
from popstar import atmospheres as atm
from scipy import interpolate
from scipy import stats
from pysynphot import spectrum
from pysynphot import ObsBandpass
from pysynphot import observation as obs
import pysynphot
from astropy import constants, units
from astropy.table import Table, Column
from popstar.imf import imf, multiplicity
from popstar.utils import objects
import pickle
import time, datetime
import math
import os, glob
import tempfile
import scipy
import matplotlib
import matplotlib.pyplot as plt
import time
import pdb

default_evo_model = evolution.MergedPisaEkstromParsec()
default_red_law = reddening.RedLawNishiyama09()
default_atm_func = atm.get_merged_atmosphere

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
    imf_in = imf.Kroupa_2001(multiplicity=multi)
    evo = evolution.MergedPisaEkstromParsec()
    atm_func = atm.get_merged_atmosphere

    log_age = 6.5
    AKs = 1.0
    distance = 8000.0
    cluster_mass = 100.

    if resolved:
        cluster = ResolvedCluster(log_age, AKs, distance, cluster_mass, imf_in, evo, atm_func)
    else:
        cluster = UnresolvedCluster(log_age, AKs, distance, cluster_mass, imf_in, evo, atm_func)

    # Plot the spectrum of the most massive star
    idx = cluster.mass_all.argmax()
    print 'Most massive star is {0:f} M_sun.'.format(cluster.mass_all[idx])
    #bigstar = cluster.spec_list_trim[idx]
    plt.figure(1)
    plt.clf()
    plt.plot(cluster.spec_list_trim[idx]._wavetable, cluster.spec_list_trim[idx]._fluxtable, 'k.')

    # Plot an integrated spectrum of the whole cluster.
    wave, flux = cluster.spec_trim._wavetable, cluster.spec_trim._fluxtable
    plt.figure(2)
    plt.clf()
    plt.plot(wave, flux, 'k.')

    return

class Cluster(object):
    def __init__(self, logAge, AKs, distance, clusterMass,
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
        results = imf.generate_cluster(clusterMass)
        
        self.mass = results[0] # masses of the stars in the simulated cluster
        isMultiple = results[1]
        compMasses = results[2]
        systemMasses = results[3]

        return

    
class ResolvedCluster(Cluster):
    def __init__(self, logAge, AKs, distance, clusterMass,
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
    def __init__(self, logAge, AKs, distance, clusterMass,
                 imf, evolution_model, atmosphere_func,
                 wave_range=[19000.,25000.],
                 verbose=False): 

        Cluster.__init__(self, logAge, AKs, distance, clusterMass,
                 imf, evolution_model, atmosphere_func,
                 verbose=False)
        
        Iso = Isochrone(logAge, AKs, distance, evolution_model)

        temp = np.zeros(len(self.mass), dtype=float)
        self.mass_all = np.zeros(len(self.mass), dtype=float)
        self.spec_list = [None] * len(self.mass)
        self.spec_list_trim = [None] * len(self.mass)

        t1 = time.time()
        for ii in range(len(self.mass)):
            # Find the closest model mass (returns None, if nothing with dm = 0.1
            mdx = match_model_mass(Iso.points['mass'],self.mass[ii])
            if mdx == None:
                continue

            temp[ii] = Iso.points['Teff'][mdx]
            self.mass_all[ii] = Iso.points['mass'][mdx]
            self.spec_list[ii] = Iso.spec_list[mdx]
            self.spec_list_trim[ii] = spectrum.trimSpectrum(Iso.spec_list[mdx],wave_range[0],wave_range[1])

        t2 = time.time()
        print 'Mass matching took {0:f} s.'.format(t2-t1)

        # Get rid of the bad ones
        idx = np.where(temp != 0)[0]
        cdx = np.where(temp == 0)[0]

        self.mass_all = self.mass_all[idx]
        self.spec_list = [self.spec_list[iidx] for iidx in idx]
        self.spec_list_trim = [self.spec_list_trim[iidx] for iidx in idx]
        
        self.spec_tot_full = np.sum(self.spec_list)
        
        self.spec_trim = spectrum.trimSpectrum(self.spec_tot_full,wave_range[0],wave_range[1])

        self.mass_tot = np.sum(self.mass_all)
        print 'Total mass is {0:f} M_sun'.format(self.mass_tot)
        self.log_age = logAge

        return

        
def get_evo_model_by_string(evo_model_string):
    return getattr(evolution, evo_model_string)

class Isochrone(object):
    def __init__(self, logAge, AKs, distance,
                 evo_model=default_evo_model, atm_func=default_atm_func,
                 redlaw=default_red_law, mass_sampling=2,
                 wave_range=[5000, 42500]):
        """
        Parameters
        ----------
        logAge : float
        AKs : float
        distance : float
        evModel : model cl
        mass_sampling - Sample the raw isochrone every ## steps. The default
                       is mass_sampling = 10, which takes every 10th point.
                       The isochrones are already very finely sampled. Must be
                       an integer value.
        wave_range : list
            length=2 list with the wavelength min/max of the final spectra.
            Units are Angstroms.
        """

        t1 = time.time()
        
        c = constants

        # Get solar metallicity models for a population at a specific age.
        evol = evo_model.isochrone(age=10**logAge)  # solar metallicity

        #Eliminate cases where log g is less than 0
        idx = np.where(evol['logg'] > 0)
        evol = evol[idx]

        # Trim down the table by selecting every Nth point where
        # N = mass sampling factor.
        evol = evol[::mass_sampling]

        # Determine which stars are WR stars.
        evol['isWR'] = evol['logT'] != evol['logT_WR']

        # Give luminosity, temperature, mass, radius units (astropy units).
        L_all = 10**evol['logL'] * c.L_sun # luminsoity in erg/s
        T_all = 10**evol['logT'] * units.K
        R_all = np.sqrt(L_all / (4.0 * math.pi * c.sigma_sb * T_all**4))
        mass_all = evol['mass'] * units.Msun # masses in solar masses
        logg_all = evol['logg']

        # Define the table that contains the "average" properties for each star.
        tab = Table([L_all, T_all, R_all, mass_all, logg_all], names=['L', 'Teff', 'R', 'mass', 'logg'])

        # Initialize output for stellar spectra
        self.spec_list = []

        # For each temperature extract the synthetic photometry.
        for ii in range(len(tab['Teff'])):
            gravity = float( tab['logg'][ii] )
            L = float( tab['L'].quantity[ii].cgs / (units.erg / units.s)) # in erg/s
            T = float( tab['Teff'].quantity[ii] / units.K)               # in Kelvin
            R = float( tab['R'].quantity[ii].to('pc') / units.pc)              # in pc

            # Get the atmosphere model now. Wavelength is in Angstroms
            star = atm_func(temperature=T, gravity=gravity)
            
            # Trim wavelength range down to JHKL range (0.5 - 4.25 microns)
            star = spectrum.trimSpectrum(star, wave_range[0], wave_range[1])

            # Convert into flux observed at Earth (unreddened)
            star *= (R / distance)**2  # in erg s^-1 cm^-2 A^-1
            red = redlaw.reddening(AKs).resample(star.wave)  ## check this
            star *= red
            self.spec_list.append(star)

        tab.meta['RedLaw'] = redlaw.name
        tab.meta['AtmFunc'] = atm_func.__name__
        tab.meta['EvoModel'] = type(evo_model).__name__
        tab.meta['LogAge'] = logAge
        tab.meta['AKs'] = AKs
        tab.meta['distance'] = distance
        tab.meta['wave_min'] = wave_range[0]
        tab.meta['wave_max'] = wave_range[1]

        self.points = tab

        t2 = time.time()
        print 'Isochrone generation took {0:f} s.'.format(t2-t1)
        
        return

    def trim(self, keep_indices):
        # Convert luminosity to erg/s
        self.points = self.points[keep_indices]
        self.spec_list = self.spec_list[keep_indices]

        return

class IsochronePhot(Isochrone):
    def __init__(self, logAge, AKs, distance,
                 evo_model=default_evo_model, atm_func=default_atm_func,
                 redlaw=default_red_law, mass_sampling=2, iso_dir='./',
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

        """Make an isochrone with photometry in various filters.

        Description
        -----------
        Make an isochrone with photometry in various filters. Load from file
        or save to file if possible.

        Parameters
        ----------

        Returns
        -------
                 
        """
        
        Isochrone.__init__(self, logAge, AKs, distance,
                           evo_model=evo_model, atm_func=atm_func,
                           redlaw=red_law, mass_sampling=mass_sampling)

        # Make and input/output file name for the stored isochrone photometry.
        inFileFmt = '{0}iso_{1:.2f}_{2:4.2f}_{3:4s}.fits'
        inFile = inFileFmt.format(iso_dir, logAge, AKs, str(distance).zfill(5))

        if not os.path.exists(inFile):
            make_observed_isochrone_hst(logAge=logAge, AKs=AKs, distance=distance,
                                        iso_dir=iso_dir, massSampling=massSampling,
                                        filters=filters)
        else:
            iso = Table.read(inFile)
            # Add some error checking.

        return

    def make_photometry(self):
        """
        Make synthetic photometry for the specified filters. This function
        udpates the self.points table to include new columns with the
        photometry.

        Parameters
        ----------
        filters : dictionary
            A dictionary containing the filter name (for the output columns)
            and the filter specification string that can be processed by pysynphot.

        
        """
        startTime = time.time()

        print 'Making photometry for isochrone: log(t) = %.2f  AKs = %.2f  dist = %d' % \
            (logAge, AKs, distance)
        print '     Starting at: ', datetime.datetime.now(), '  Usually takes ~5 minutes'

        npoints = len(self.points)
        verbose_fmt = 'M = {0:7.3f} Msun  T = {1:5d} K  m_{2:s} = {3:4.2f}'

        # Loop through the filters, get filter info, make photometry for
        # all stars in this filter.
        for filt_name, filt_str in self.filters.iteritems():
            prt_fmt = 'Starting filter: {0:s}   Elapse time: {1:d} seconds'
            print prt_fmt.format(filt_name, time.time() - startTime)
            
            filt = get_filter_info(filt_str)

            # Make the column to hold magnitudes in this filter. Add to points table.
            col_name = 'mag' + filt_name
            mag_col = Column(np.zeros(npoints, dtype=float), name=col_name)
            self.points.add_column(mag_col)
            
            # Loop through each star in the isochrone and do the filter integration
            for ss in range(npoints):
                star = self.spec_list[ss]  # These are already extincted, observed spectra.
                star_mag = mag_in_filter(star, filt)
                
                self.points[col_name][ss] = star_mag
        
                if (self.verbose and (ss % 100) == 0):
                    print verbose_mt.format(self.points['mass'][ss], self.points['Teff'][ss],
                                            filtName, star_mag)

        endTime = time.time()
        print '      Time taken: %d seconds' % (endTime - startTime)

        return

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
        from nirc2 import synthetic as nirc2syn
    
        tmp = name.split(',')
        filterName = tmp[-1]
        filt = nirc2syn.FilterNIRC2(filterName)
    else:
        filt = ObsBandpass(name)

    vega_obs = obs.Observation(vega, filt, binset=filt.wave, force='taper')
    vega_flux = vega_obs.binflux.sum()
    vega_mag = 0.03

    filt.flux0 = vega_flux
    filt.mag0 = vega_mag
    
    return filt

# Little helper utility to get the magnitude of an object through a filter.
def mag_in_filter(star, filt):
    """
    Assumes that extinction is already resampled to same wavelengths
    as filter.
    """
    star_in_filter = obs.Observation(star, filt, binset=filt.wave, force='taper')
    star_flux = star_in_filter.binflux.sum()
    star_mag = -2.5 * math.log10(star_flux / filt.flux0) + filt.mag0

    return star_mag

# def model_young_cluster_new(logAge, AKs, distance, iso_dir,
#                             evModel='mergedPisaEkstromParsec',
#                             imfSlopes=np.array([-2.35]), massLimits=np.array([1, 150]),
#                             clusterMass=1e4, makeMultiples=False,
#                             MFamp=defaultMFamp, MFindex=defaultMFindex,
#                             CSFamp=defaultCSFamp, CSFindex=defaultCSFindex,
#                             CSFmax=defaultCSFmax,
#                             qMin=0.01, qIndex=-0.4, verbose=False):
#     """
#     Code to model a cluster with user-specified logAge, AKs, and distance.
#     Must also specify directory containing the isochrone (made using popstar
#     synthetic code).

#     Can also specify IMF slope, mass limits, cluster mass, and parameters for
#     multiple stars
#     """
#     c = constants

#     logAgeString = '0%d' % (int(logAge * 100))

#     iso = load_isochrone(logAge=logAge, AKs=AKs, distance=distance, evModel=evModel,
#                                    iso_dir=iso_dir)

#     # Sample a power-law IMF randomly
#     results = imf.sample_imf(massLimits, imfSlopes, clusterMass,
#                              makeMultiples=makeMultiples,
#                              multiMFamp=MFamp, multiMFindex=MFindex,
#                              multiCSFamp=CSFamp, multiCSFindex=CSFindex,
#                              multiCSFmax=CSFmax,
#                              multiQmin=qMin, multiQindex=qIndex)
#     mass = results[0]
#     isMultiple = results[1]
#     compMasses = results[2]
#     systemMasses = results[3]

#     mag814w = np.zeros(len(mass), dtype=float)
#     mag127m = np.zeros(len(mass), dtype=float)
#     mag139m = np.zeros(len(mass), dtype=float)
#     mag153m = np.zeros(len(mass), dtype=float)
#     magJ = np.zeros(len(mass), dtype=float)
#     magH = np.zeros(len(mass), dtype=float)
#     magK = np.zeros(len(mass), dtype=float)
#     magKp = np.zeros(len(mass), dtype=float)
#     magL = np.zeros(len(mass), dtype=float)
#     temp = np.zeros(len(mass), dtype=float)
#     logg = np.zeros(len(mass), dtype=float)
#     logL = np.zeros(len(mass), dtype=float)
#     isWR = np.zeros(len(mass), dtype=bool)
    
#     def match_model_mass(theMass):
#         dm = np.abs(iso['mass'] - theMass)
#         mdx = dm.argmin()

#         # Model mass has to be within 2% of the desired mass
#         if (dm[mdx] / theMass) > 0.1:
#             return None
#         else:
#             return mdx
        
#     for ii in range(len(mass)):
#         # Find the closest model mass (returns None, if nothing with dm = 0.1
#         mdx = match_model_mass(mass[ii])
#         if mdx == None:
#             continue

#         mag814w[ii] = iso['mag814w'][mdx]
#         mag127m[ii] = iso['mag127m'][mdx]
#         mag139m[ii] = iso['mag139m'][mdx]
#         mag153m[ii] = iso['mag153m'][mdx]
#         magJ[ii] = iso['magJ'][mdx]
#         magH[ii] = iso['magH'][mdx]
#         magK[ii] = iso['magK'][mdx]
#         magKp[ii] = iso['magKp'][mdx]
#         magL[ii] = iso['magL'][mdx]
        
#         temp[ii] = iso['logT'][mdx]
#         logg[ii] = iso['logg'][mdx]
#         logL[ii] = iso['logL'][mdx]
#         isWR[ii] = iso['isWR'][mdx]


#         # Determine if this system is a binary.
#         if isMultiple[ii]:
#             n_stars = len(compMasses[ii])
#             for cc in range(n_stars):
#                 mdx_cc = match_model_mass(compMasses[ii][cc])
#                 if mdx_cc != None:
#                     f1 = 10**(-mag[ii]/2.5)
#                     f2 = 10**(-iso.mag[mdx_cc]/2.5)
#                     mag[ii] = -2.5 * np.log10(f1 + f2)
#                 else:
#                     print 'Rejected a companion %.2f' % compMasses[ii][cc]
        

#     # Get rid of the bad ones
#     idx = np.where(temp != 0)[0]
#     cdx = np.where(temp == 0)[0]

#     if len(cdx) > 0 and verbose:
#         print 'Found %d stars out of mass range: Minimum bad mass = %.1f' % \
#             (len(cdx), mass[cdx].min())

#     mass = mass[idx]
#     mag814w = mag814w[idx]
#     mag127m = mag127m[idx]
#     mag139m = mag139m[idx]
#     mag153m = mag153m[idx]
#     magJ = magJ[idx]
#     magH = magH[idx]
#     magK = magK[idx]
#     magKp = magKp[idx]
#     magL = magL[idx]
    
    
#     temp = temp[idx]
#     logg = logg[idx]
#     logL = logL[idx]
#     isWR = isWR[idx]
    
#     isMultiple = isMultiple[idx]
#     systemMasses = systemMasses[idx]
#     if makeMultiples:
#         compMasses = [compMasses[ii] for ii in idx]

#     idx_noWR = np.where(isWR == False)[0]

#     mag127m_noWR = mag127m[idx_noWR]
#     num_WR = len(mag127m) - len(idx_noWR)

#     cluster = objects.DataHolder()
#     cluster.mass = mass
#     cluster.Teff = temp
#     cluster.logg = logg
#     cluster.logL = logL
#     cluster.isWR = isWR
#     cluster.mag814w = mag814w
#     cluster.mag127m = mag127m
#     cluster.mag139m = mag139m
#     cluster.mag153m = mag153m
#     cluster.magJ = magJ
#     cluster.magH = magH
#     cluster.magK = magK
#     cluster.magKp = magKp
#     cluster.magL = magL
#     cluster.isMultiple = isMultiple
#     cluster.compMasses = compMasses
#     cluster.systemMasses = systemMasses

#     cluster.idx_noWR = idx_noWR
#     cluster.mag127m_noWR = mag127m_noWR
#     cluster.num_WR = num_WR

#     # Summary parameters
#     cluster.logAge = logAge
#     cluster.AKs = AKs
#     cluster.distance = distance
#     cluster.massLimits = massLimits
#     cluster.imfSlopes = imfSlopes
#     cluster.sumIMFmass = clusterMass
#     cluster.makeMultiples = makeMultiples
#     cluster.MFamp = MFamp
#     cluster.MFindex = MFindex
#     cluster.CSFamp = CSFamp
#     cluster.CSFindex = CSFindex
#     cluster.CSFmax = CSFmax
#     cluster.qMin = qMin
#     cluster.qIndex = qIndex
            
#     return cluster, iso
