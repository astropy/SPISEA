import numpy as np
import pylab as plt
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

class Cluster(object):
    def __init__(self, iso, imf, cluster_mass, verbose=False):
        """
        Code to model a cluster with user-specified logAge, AKs, and distance.
        Must also specify directory containing the isochrone (made using popstar
        synthetic code).

        Can also specify IMF slope, mass limits, cluster mass, and parameters for
        multiple stars
        """
        self.verbose = verbose
        self.iso = iso
        self.imf = imf
        self.cluster_mass = cluster_mass

        return

    
class ResolvedCluster(Cluster):
    def __init__(self, iso, imf, cluster_mass, verbose=False):
        # Doesn't do much.
        Cluster.__init__(self, iso, imf, cluster_mass, verbose=verbose)
        
        # Sample the IMF to build up our cluster mass.
        mass, isMulti, compMass, sysMass = imf.generate_cluster(cluster_mass)

        ##### 
        # Make a table to contain all the information about each stellar system.
        #####
        star_systems = Table([mass, isMulti, sysMass],
                             names=['mass', 'isMultiple', 'systemMass'])
        N_systems = len(star_systems)

        # Add columns for the Teff, L, logg, isWR for the primary stars.
        star_systems.add_column( Column(np.zeros(N_systems, dtype=float), name='Teff') )
        star_systems.add_column( Column(np.empty(N_systems, dtype=float), name='L') )
        star_systems.add_column( Column(np.empty(N_systems, dtype=float), name='logg') )
        star_systems.add_column( Column(np.empty(N_systems, dtype=float), name='isWR') )

        # Add the filter columns to the table. They are empty so far.
        # Keep track of the filter names in : filt_names
        filt_names = []
        for col_name in iso.points.colnames:
            if 'mag' in col_name:
                filt_names.append(col_name)
                mag_col = Column(np.empty(N_systems, dtype=float), name=col_name)
                star_systems.add_column(mag_col)

        # Loop through each of the systems and assign magnitudes.
        for ii in range(N_systems):
            # Find the closest model mass (returns None, if nothing with dm = 0.1
            mdx = match_model_mass(iso.points['mass'], star_systems['mass'][ii])
            if mdx == None:
                # print 'Rejected a primary star {0:.2f}'.format(star_systems['mass'][ii])
                continue

            star_systems['Teff'][ii] = iso.points['Teff'][mdx]
            star_systems['L'][ii] = iso.points['L'][mdx]
            star_systems['logg'][ii] = iso.points['logg'][mdx]
            star_systems['isWR'][ii] = iso.points['isWR'][mdx]
            
            for filt in filt_names:
                star_systems[filt][ii] = iso.points[filt][mdx]
                
        # Get rid of the bad ones
        idx = np.where(star_systems['Teff'] > 0)[0]
        if len(idx) != N_systems and verbose:
            print 'Found {0:d} stars out of mass range'.format(N_systems - len(idx))

        star_systems = star_systems[idx]
        N_systems = len(star_systems)

        #####
        #    MULTIPLICITY                 
        # Make a second table containing all the companion-star masses.
        # This table will be much longer... here are the arrays:
        #    sysIndex - the index of the system this star belongs too
        #    mass - the mass of this individual star.
        if imf.make_multiples:
            # Clean up companion stuff (which we haven't handeled yet)
            compMass = [compMass[ii] for ii in idx]
            N_companions = np.array([len(star_masses) for star_masses in compMass])
            star_systems.add_column( Column(N_companions, name='N_companions') )

            N_comp_tot = N_companions.sum()
            system_index = np.repeat(np.arange(N_systems), N_companions)

            companions = Table([system_index], names=['system_idx'])

            # Add columns for the Teff, L, logg, isWR, filters for the companion stars.
            companions.add_column( Column(np.zeros(N_comp_tot, dtype=float), name='mass') )
            companions.add_column( Column(np.zeros(N_comp_tot, dtype=float), name='Teff') )
            companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name='L') )
            companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name='logg') )
            companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name='isWR') )
            for filt in filt_names:
                companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name=filt) )
            
            kk = 0

            # Loop through each star system
            for ii in range(N_systems):
                # Determine if this system is a multiple star system.
                if star_systems['isMultiple'][ii]:
                    
                    # Loop through the companions in this system
                    for cc in range(N_companions[ii]):
                        companions['mass'][kk] = compMass[ii][cc]
                        mdx_cc = match_model_mass(iso.points['mass'], compMass[ii][cc])
                        
                        if mdx_cc != None:
                            companions['Teff'][kk] = iso.points['Teff'][mdx_cc]
                            companions['L'][kk] = iso.points['L'][mdx_cc]
                            companions['logg'][kk] = iso.points['logg'][mdx_cc]
                            companions['isWR'][kk] = iso.points['isWR'][mdx_cc]
                            
                            for filt in filt_names:
                                f1 = 10**(-star_systems[filt][ii] / 2.5)
                                f2 = 10**(-iso.points[filt][mdx_cc] / 2.5)
                                
                                companions[filt][kk] = iso.points[filt][mdx_cc]
                                star_systems[filt][ii] = -2.5 * np.log10(f1 + f2)

                            kk += 1
                        
                        # else:
                        #     print 'Rejected a companion %.2f' % compMass[ii][cc]
                            
            # Get rid of the bad ones
            idx = np.where(companions['Teff'] > 0)[0]
            if len(idx) != N_comp_tot and self.verbose:
                print 'Found {0:d} companions out of mass range'.format(N_comp_tot - len(idx))

            # Double check that everything behaved properly.
            assert companions['mass'][idx].min() > 0

        # Save our arrays to the object
        self.star_systems = star_systems
        
        if imf.make_multiples:
            self.companions = companions

        return

            
class UnresolvedCluster(Cluster):
    def __init__(self, iso, imf, cluster_mass,
                 wave_range=[5000, 50000], verbose=False):
        # Doesn't do much.
        Cluster.__init__(self, iso, imf, cluster_mass, verbose=verbose)
        
        # Sample a power-law IMF randomly
        mass, isMulti, compMass, sysMass = imf.generate_cluster(cluster_mass)
        
        temp = np.zeros(len(mass), dtype=float)
        spec_list = [None] * len(mass)

        t1 = time.time()
        for ii in range(len(mass)):
            # Find the closest model mass (returns None, if nothing with dm = 0.1
            mdx = match_model_mass(iso.points['mass'], mass[ii])
            if mdx == None:
                continue

            temp[ii] = iso.points['Teff'][mdx]
            spec_list[ii] = iso.spec_list[mdx]

        t2 = time.time()
        print 'Mass matching took {0:f}s'.format(t2-t1)

        # Get rid of the bad ones
        idx = np.where(temp != 0)[0]
        cdx = np.where(temp == 0)[0]

        spec_list = [spec_list[iidx] for iidx in idx]
        
        self.spec_tot_full = np.sum(spec_list, 0)
        print 'Spec summing took {0:f}s'.format(t2-t1)
        
        self.spec_trim = spectrum.trimSpectrum(self.spec_tot_full,
                                               wave_range[0], wave_range[1])
        print 'Spec trimming took {0:f}s'.format(t2-t1)

        self.mass_all = np.sum(sysMass[idx])

        return

        
def get_evo_model_by_string(evo_model_string):
    return getattr(evolution, evo_model_string)

class Isochrone(object):
    def __init__(self, logAge, AKs, distance,
                 evo_model=default_evo_model, atm_func=default_atm_func,
                 red_law=default_red_law, mass_sampling=1,
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

        c = constants

        # Get solar metallicity models for a population at a specific age.
        # Takes about 0.1 seconds.
        evol = evo_model.isochrone(age=10**logAge)  # solar metallicity

        # Eliminate cases where log g is less than 0
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
        isWR_all = evol['isWR']

        # Define the table that contains the "average" properties for each star.
        tab = Table([L_all, T_all, R_all, mass_all, logg_all, isWR_all],
                    names=['L', 'Teff', 'R', 'mass', 'logg', 'isWR'])

        # Initialize output for stellar spectra
        self.spec_list = []

        # For each temperature extract the synthetic photometry.
        for ii in range(len(tab['Teff'])):
            # Loop is currently taking about 0.11 s per iteration

            gravity = float( logg_all[ii] )
            L = float( L_all[ii].cgs / (units.erg / units.s)) # in erg/s
            T = float( T_all[ii] / units.K)               # in Kelvin
            R = float( R_all[ii].to('pc') / units.pc)              # in pc

            # Get the atmosphere model now. Wavelength is in Angstroms
            # This is the time-intensive call... everything else is negligable.
            star = atm_func(temperature=T, gravity=gravity)
            
            # Trim wavelength range down to JHKL range (0.5 - 4.25 microns)
            star = spectrum.trimSpectrum(star, wave_range[0], wave_range[1])

            # Convert into flux observed at Earth (unreddened)
            star *= (R / distance)**2  # in erg s^-1 cm^-2 A^-1

            # Redden the spectrum. This doesn't take much time at all.
            red = red_law.reddening(AKs).resample(star.wave) 
            star *= red

            # Save the final spectrum to our spec_list for later use.            
            self.spec_list.append(star)

        # Append all the meta data to the summary table.
        tab.meta['REDLAW'] = red_law.name
        tab.meta['ATMFUNC'] = atm_func.__name__
        tab.meta['EVOMODEL'] = type(evo_model).__name__
        tab.meta['LOGAGE'] = logAge
        tab.meta['AKS'] = AKs
        tab.meta['DISTANCE'] = distance
        tab.meta['WAVEMIN'] = wave_range[0]
        tab.meta['WAVEMAX'] = wave_range[1]

        self.points = tab

        return

    def trim(self, keep_indices):
        # Convert luminosity to erg/s
        self.points = self.points[keep_indices]
        self.spec_list = self.spec_list[keep_indices]

        return

    def plot_HR_diagram(self, savefile=None):
        """
        Make a standard HR diagram for this isochrone.
        """
        plt.clf()
        plt.loglog(self.points['Teff'], self.points['L'],
                   color='black', linestyle='solid', marker='+')
        plt.gca().invert_xaxis()
        plt.xlabel(r'T$_{\mathrm{eff}}$ (K)')
        plt.ylabel('Luminosity (erg / s)')
        
        fmt_title = 'logAge={0:.2f}, d={1:.2f} kpc, AKs={2:.2f}'
        plt.title(fmt_title.format(self.points.meta['LOGAGE'],
                                  self.points.meta['DISTANCE']/1e3,
                                  self.points.meta['AKS']))

        if savefile != None:
            plt.savefig(savefile)
        
        return

    def plot_mass_luminosity(self, savefile=None):
        """
        Make a standard mass-luminosity relation plot for this isochrone.
        """
        plt.clf()
        plt.loglog(self.points['mass'], self.points['L'], 'k.')
        plt.xlabel(r'Mass (M$_\odot$)')
        plt.ylabel('Luminosity (erg / s)')
        
        fmt_title = 'logAge={0:.2f}, d={1:.2f} kpc, AKs={2:.2f}'
        plt.title(fmt_title.format(self.points.meta['LOGAGE'],
                                  self.points.meta['DISTANCE']/1e3,
                                  self.points.meta['AKS']))
        
        if savefile != None:
            plt.savefig(savefile)
        
        return

class IsochronePhot(Isochrone):
    def __init__(self, logAge, AKs, distance,
                 evo_model=default_evo_model, atm_func=default_atm_func,
                 red_law=default_red_law, mass_sampling=1, iso_dir='./',
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
        
        # Make and input/output file name for the stored isochrone photometry.
        save_file_fmt = '{0}iso_{1:.2f}_{2:4.2f}_{3:4s}.fits'
        self.save_file = save_file_fmt.format(iso_dir, logAge, AKs, str(distance).zfill(5))

        # Expected filters
        self.filters = filters

        if not os.path.exists(self.save_file):
            Isochrone.__init__(self, logAge, AKs, distance,
                               evo_model=evo_model, atm_func=atm_func,
                               red_law=red_law, mass_sampling=mass_sampling)
            self.verbose = True
            self.make_photometry()
        else:
            self.points = Table.read(self.save_file)
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

        meta = self.points.meta

        print 'Making photometry for isochrone: log(t) = %.2f  AKs = %.2f  dist = %d' % \
            (meta['LOGAGE'], meta['AKS'], meta['DISTANCE'])
        print '     Starting at: ', datetime.datetime.now(), '  Usually takes ~5 minutes'

        npoints = len(self.points)
        verbose_fmt = 'M = {0:7.3f} Msun  T = {1:5.0f} K  m_{2:s} = {3:4.2f}'

        # Loop through the filters, get filter info, make photometry for
        # all stars in this filter.
        for filt_name, filt_str in self.filters.iteritems():
            prt_fmt = 'Starting filter: {0:s}   Elapsed time: {1:.2f} seconds'
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
                    print verbose_fmt.format(self.points['mass'][ss], self.points['Teff'][ss],
                                             filt_name, star_mag)

        endTime = time.time()
        print '      Time taken: {0:.2f} seconds'.format(endTime - startTime)

        if self.save_file != None:
            self.points.write(self.save_file)

        return

    def plot_CMD(self, mag1, mag2, savefile=None):
        """
        Make a CMD with mag1 vs. mag1 - mag2

        Parameters
        ----------
        mag1 : string
            The name of the first magnitude column to be plotted.
        mag2 : string
            The name of the second magnitude column to be plotted.
        savefile : string (default None)
            If a savefile is specified, then the plot will be saved to that file. 
        """
        plt.clf()
        plt.plot(self.points[mag1] - self.points[mag2], self.points[mag1],
                 color='black', linestyle='solid', marker='+')
        plt.gca().invert_yaxis()
        plt.xlabel(mag1 + ' - ' + mag2 + ' (mag)')
        plt.ylabel(mag1 + ' (mag)')
        
        fmt_title = 'logAge={0:.2f}, d={1:.2f} kpc, AKs={2:.2f}'
        plt.title(fmt_title.format(self.points.meta['LOGAGE'],
                                  self.points.meta['DISTANCE']/1e3,
                                  self.points.meta['AKS']))

        if savefile != None:
            plt.savefig(savefile)
        
        return

    def plot_mass_magnitude(self, mag, savefile=None):
        """
        Make a standard mass-luminosity relation plot for this isochrone.
        
        Parameters
        ----------
        mag : string
            The name of the magnitude column to be plotted.
        savefile : string (default None)
            If a savefile is specified, then the plot will be saved to that file. 
        """
        plt.clf()
        plt.semilogx(self.points['mass'], self.points[mag], 'k.')
        plt.gca().invert_yaxis()
        plt.xlabel(r'Mass (M$_\odot$)')
        plt.ylabel(mag + ' (mag)')
        
        fmt_title = 'logAge={0:.2f}, d={1:.2f} kpc, AKs={2:.2f}'
        plt.title(fmt_title.format(self.points.meta['LOGAGE'],
                                  self.points.meta['DISTANCE']/1e3,
                                  self.points.meta['AKS']))
        
        if savefile != None:
            plt.savefig(savefile)
        
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
        
        # Convert to ArraySpectralElement for resampling.
        filt = spectrum.ArraySpectralElement(filt.wave, filt.throughput,
                                             waveunits=filt.waveunits,
                                             name=filt.name)

    # Resample the filter to have 1500 points across. More is excessive.
    if len(filt.wave) > 1500:
        idx = np.where(filt.throughput > 0.001)[0]
        new_wave = np.linspace(filt.wave[idx[0]], filt.wave[idx[-1]], 1500, dtype=float)
        filt = filt.resample(new_wave) 

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
