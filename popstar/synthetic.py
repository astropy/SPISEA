import numpy as np
import pylab as plt
from popstar import reddening
from popstar import evolution
from popstar import atmospheres as atm
from popstar import filters 
from scipy import interpolate
from scipy import stats
from scipy.special import erf
from pysynphot import spectrum
from pysynphot import ObsBandpass
from pysynphot import observation as obs
import pysynphot
from astropy import constants, units
from astropy.table import Table, Column, MaskedColumn
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
import warnings
import pdb
from scipy.spatial import cKDTree as KDTree

default_evo_model = evolution.MergedBaraffePisaEkstromParsec()
default_red_law = reddening.RedLawNishiyama09()
default_atm_func = atm.get_merged_atmosphere

def Vega():
    # Use Vega as our zeropoint... assume V=0.03 mag and all colors = 0.0
    # These parameters are defined in Girardi+02
    vega = atm.get_kurucz_atmosphere(temperature=9550, 
                                     gravity=3.95,
                                     metallicity=-0.5)

    vega = spectrum.trimSpectrum(vega, 3000, 50000)

    # This is (R/d)**2 as reported by Girardi et al. 2002, page 198, col 1.
    # and is used to convert to flux observed at Earth.
    vega *= 6.247e-17 
    
    return vega

vega = Vega()

class Cluster(object):
    def __init__(self, iso, imf, cluster_mass, ifmf=None, verbose=False):
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
        self.ifmf = ifmf
        self.cluster_mass = cluster_mass

        return

class ResolvedCluster(Cluster):
    def __init__(self, iso, imf, cluster_mass, ifmf=None, save_dir='./', verbose=True):
        # Save to object variables
        Cluster.__init__(self, iso, imf, cluster_mass, ifmf=ifmf, verbose=verbose)

        t1 = time.time()
        ##### 
        # Sample the IMF to build up our cluster mass.
        #####
        mass, isMulti, compMass, sysMass = imf.generate_cluster(cluster_mass)

        # Figure out the filters we will make.
        self.filt_names = self.set_filter_names()

        ##### 
        # Make a table to contain all the information about each stellar system.
        #####
        #t3 = time.time()
        star_systems = self._make_star_systems_table_interp(mass, isMulti, sysMass)
        #t4 = time.time()
        #print 'make star systems: {0}'.format(t4 - t3)
        
        # Trim out bad systems; specifically, stars with masses outside those provided
        # by the model isochrone
        star_systems, compMass = self._remove_bad_systems(star_systems, compMass, self.ifmf)

        ###
        # Calculate remnant masses and identity using ifmf, if desired
        ###
        if self.ifmf != None:
            # Initialize remnant ID, remnant_mass columns. Default ID is 0, mass = -99 for
            # non-compact objects
            remnant_id = np.zeros(len(star_systems['mass']))
            remnant_mass = np.ones(len(star_systems['mass'])) * -99

            # Identify compact objects as those with Teff = 0 (indicating they are not
            # present in isochrone)
            idx_rem = np.where(star_systems['Teff'] == 0)
            
            # Calculate remnant mass and ID for compact objects; update remnant_id and
            # remnant_mass arrays accordingly
            r_mass_tmp, r_id_tmp = self.ifmf.generate_death_mass_distribution(star_systems['mass'][idx_rem])
            remnant_mass[idx_rem] = r_mass_tmp
            remnant_id[idx_rem] = r_id_tmp

            # Mask remnant mass where it is not relevant (e.g. not a compact object or
            # outside mass range IFMF is defined for)
            remnant_mass = np.ma.masked_where((remnant_id <= 0), remnant_mass)
            
            # Add columns to star_systems table
            remnant_mass_col = MaskedColumn(remnant_mass, name='Rem_mass')
            remnant_id_col = Column(remnant_id, name='Rem_ID')
            star_systems.add_column(remnant_id_col)
            star_systems.add_column(remnant_mass_col)
            
        ##### 
        # Make a table to contain all the information about companions.
        #####
        if self.imf.make_multiples:
            #t5 = time.time()
            companions = self._make_companions_table_interp(star_systems, compMass, self.ifmf)
            #t6 = time.time()
            #print 'make comp systems: {0}'.format(t6 - t5)
        #####
        # Save our arrays to the object
        #####
        self.star_systems = star_systems
        
        if self.imf.make_multiples:
            self.companions = companions

        return

    def set_filter_names(self):
        filt_names = []
        
        for col_name in self.iso.points.colnames:
            if 'm_' in col_name:
                filt_names.append(col_name)

        return filt_names
        

    def _make_star_systems_table(self, mass, isMulti, sysMass):
        """
        Make a star_systems table and get synthetic photometry for each primary star.
        """
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
        for filt in self.filt_names:
            star_systems.add_column( Column(np.empty(N_systems, dtype=float), name=filt) )

        mdx_all = match_model_masses(self.iso.points['mass'], star_systems['mass'])
        idx_good = np.where(mdx_all >= 0)[0]
        mdx_good = mdx_all[idx_good]

        if len(idx_good) != len(mdx_all):
            msg = 'Rejected {0:d} of {1:d} primary stars' 
            print( msg.format(len(mdx_all) - len(idx_good), N_systems))
            foo = np.where(mdx_all == -1)[0]
            print( star_systems['mass'][foo])

        
        star_systems['Teff'][idx_good] = self.iso.points['Teff'][mdx_good]
        star_systems['L'][idx_good] = self.iso.points['L'][mdx_good]
        star_systems['logg'][idx_good] = self.iso.points['logg'][mdx_good]
        star_systems['isWR'][idx_good] = self.iso.points['isWR'][mdx_good]

        for filt in self.filt_names:
            star_systems[filt][idx_good] = self.iso.points[filt][mdx_good]
    
        return star_systems

    def _make_star_systems_table_interp(self, mass, isMulti, sysMass):
        """
        Make a star_systems table and get synthetic photometry for each primary star.
        """
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
        for filt in self.filt_names:
            star_systems.add_column( Column(np.empty(N_systems, dtype=float), name=filt) )

        iso_pts = self.iso.points
        interp_Teff = interpolate.interp1d(iso_pts['mass'], iso_pts['Teff'], kind='linear', bounds_error=False, fill_value=0)
        interp_L = interpolate.interp1d(iso_pts['mass'], iso_pts['L'], kind='linear', bounds_error=False, fill_value=0)
        interp_logg = interpolate.interp1d(iso_pts['mass'], iso_pts['logg'], kind='linear', bounds_error=False, fill_value=0)
        interp_isWR = interpolate.interp1d(iso_pts['mass'], iso_pts['isWR'], kind='linear', bounds_error=False, fill_value=0)

        star_systems['Teff'] = interp_Teff(star_systems['mass'])
        star_systems['L'] = interp_L(star_systems['mass'])
        star_systems['logg'] = interp_logg(star_systems['mass'])
        star_systems['isWR'] = np.round(interp_isWR(star_systems['mass']))

        for filt in self.filt_names:
            interp_filt = interpolate.interp1d(iso_pts['mass'], iso_pts[filt],
                                               kind='linear',
                                               bounds_error=False, fill_value=0)
            star_systems[filt] = interp_filt(star_systems['mass'])

        return star_systems

        
    def _make_companions_table(self, star_systems, compMass):

        N_systems = len(star_systems)
        
        #####
        #    MULTIPLICITY                 
        # Make a second table containing all the companion-star masses.
        # This table will be much longer... here are the arrays:
        #    sysIndex - the index of the system this star belongs too
        #    mass - the mass of this individual star.
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
        for filt in self.filt_names:
            companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name=filt) )

        kk = 0

        # Loop through each star system
        for ii in range(N_systems):
            # Determine if this system is a multiple star system.
            if star_systems['isMultiple'][ii]:

                # Loop through the companions in this system
                for cc in range(N_companions[ii]):
                    companions['mass'][kk] = compMass[ii][cc]
                    mdx_cc = match_model_mass(self.iso.points['mass'], compMass[ii][cc])

                    if mdx_cc != None:
                        companions['Teff'][kk] = self.iso.points['Teff'][mdx_cc]
                        companions['L'][kk] = self.iso.points['L'][mdx_cc]
                        companions['logg'][kk] = self.iso.points['logg'][mdx_cc]
                        companions['isWR'][kk] = self.iso.points['isWR'][mdx_cc]

                        for filt in self.filt_names:
                            f1 = 10**(-star_systems[filt][ii] / 2.5)
                            f2 = 10**(-self.iso.points[filt][mdx_cc] / 2.5)

                            companions[filt][kk] = self.iso.points[filt][mdx_cc]
                            star_systems[filt][ii] = -2.5 * np.log10(f1 + f2)

                        kk += 1

        # Notify if we have a lot of bad ones.
        idx = np.where(companions['Teff'] > 0)[0]
        if len(idx) != N_comp_tot and self.verbose:
            print( 'Found {0:d} companions out of mass range'.format(N_comp_tot - len(idx)))

        # Double check that everything behaved properly.
        assert companions['mass'][idx].min() > 0

        return companions

    def _make_companions_table_interp(self, star_systems, compMass):

        N_systems = len(star_systems)
        
        #####
        #    MULTIPLICITY                 
        # Make a second table containing all the companion-star masses.
        # This table will be much longer... here are the arrays:
        #    sysIndex - the index of the system this star belongs too
        #    mass - the mass of this individual star.
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
        for filt in self.filt_names:
            companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name=filt) )

        # Make an array that maps system index (ii), companion index (cc) to
        # the place in the 1D companions array.
        N_comp_max = N_companions.max()
        
        comp_index = np.zeros((N_systems, N_comp_max), dtype=int)
        kk = 0
        for ii in range(N_systems):
            for cc in range(N_companions[ii]):
                comp_index[ii][cc] = kk
                kk += 1

        # Setup interpolaters.
        iso_pts = self.iso.points
        interp_Teff = interpolate.interp1d(iso_pts['mass'], iso_pts['Teff'], kind='linear', bounds_error=False, fill_value=0)
        interp_L = interpolate.interp1d(iso_pts['mass'], iso_pts['L'], kind='linear', bounds_error=False, fill_value=0)
        interp_logg = interpolate.interp1d(iso_pts['mass'], iso_pts['logg'], kind='linear', bounds_error=False, fill_value=0)
        interp_isWR = interpolate.interp1d(iso_pts['mass'], iso_pts['isWR'], kind='linear', bounds_error=False, fill_value=0)
        interp_filt = {}
        for filt in self.filt_names:
            interp_filt[filt] = interpolate.interp1d(iso_pts['mass'], iso_pts[filt],
                                                    kind='linear',
                                                    bounds_error=False, fill_value=0)
                
        # Find all the systems with at least one companion... add the flux
        # of that companion to the primary. Repeat for 2 companions,
        # 3 companions, etc.
        for cc in range(1, N_comp_max+1):
            # All systems with at least cc companions.
            idx = np.where(N_companions >= cc)[0]

            # Get the location in the companions array for each system and
            # the cc'th companion. 
            cdx = comp_index[idx, cc-1]
            
            companions['mass'][cdx] = [compMass[ii][cc-1] for ii in idx]
            comp_mass = companions['mass'][cdx]

            if len(idx) > 0:
                companions['Teff'][cdx] = interp_Teff(comp_mass)
                companions['L'][cdx] = interp_L(comp_mass)
                companions['logg'][cdx] = interp_logg(comp_mass)
                companions['isWR'][cdx] = np.round(interp_isWR(comp_mass))

                for filt in self.filt_names:
                    # Magnitude of companion
                    companions[filt][cdx] = interp_filt[filt](comp_mass)

                    # Add companion flux to system flux.
                    f1 = 10**(-star_systems[filt][idx] / 2.5)
                    f2 = 10**(-companions[filt][cdx] / 2.5)
                    star_systems[filt][idx] = -2.5 * np.log10(f1 + f2)
                    

        # Notify if we have a lot of bad ones.
        idx = np.where(companions['Teff'] > 0)[0]
        if len(idx) != N_comp_tot and self.verbose:
            print( 'Found {0:d} companions out of mass range'.format(N_comp_tot - len(idx)))

        # Double check that everything behaved properly.
        assert companions['mass'][idx].min() > 0

        return companions

    
    def _remove_bad_systems(self, star_systems, compMass, ifmf):
        """
        Helper function to remove stars with masses outside the isochrone
        mass range from the cluster. These stars are identified by having 
        a Teff = 0, as set up by _make_companions_table_interp.

        If ifmf == None, then both high and low-mass bad systems are 
        removed. If ifmf != None, then we will save the high mass systems 
        since they will be pluggedd into an ifmf later.
        """
        N_systems = len(star_systems)

        # Get rid of the bad ones
        if ifmf == None:
            idx = np.where(star_systems['Teff'] > 0)[0]
        else:
            highest_mass_iso = np.max(star_systems['mass'][np.where(star_systems['Teff'] > 0)])
            idx = np.where( (star_systems['Teff'] > 0) | (star_systems['mass'] > highest_mass_iso))[0]

        if len(idx) != N_systems and self.verbose:
            print( 'Found {0:d} stars out of mass range'.format(N_systems - len(idx)))

        star_systems = star_systems[idx]
        N_systems = len(star_systems)

        if self.imf.make_multiples:
            # Clean up companion stuff (which we haven't handeled yet)
            compMass = [compMass[ii] for ii in idx]
        
        return star_systems, compMass

class ResolvedClusterDiffRedden(ResolvedCluster):
    def __init__(self, iso, imf, cluster_mass, deltaAKs,
                 ifmf=None, red_law=default_red_law, verbose=False):

        ResolvedCluster.__init__(self, iso, imf, cluster_mass, ifmf=ifmf, verbose=verbose)

        # For a given delta_AKs (Gaussian sigma of reddening distribution at Ks),
        # figure out the equivalent delta_filt values for all other filters.
        #t1 = time.time()
        delta_red_filt = {}
        AKs = iso.points.meta['AKS']
        red_vega_lo = vega * red_law.reddening(AKs).resample(vega.wave)
        red_vega_hi = vega * red_law.reddening(AKs + deltaAKs).resample(vega.wave)

        for filt in self.filt_names:
            obs_str = get_obs_str(filt)
            filt_info = get_filter_info(obs_str)
            
            mag_lo = mag_in_filter(red_vega_lo, filt_info)
            mag_hi = mag_in_filter(red_vega_hi, filt_info)
            delta_red_filt[filt] = mag_hi - mag_lo

        # Perturb all of star systems' photometry by a random amount corresponding to
        # differential de-reddening. The distribution is normal with a width of
        # Aks +/- deltaAKs in each filter.
        rand_red = np.random.randn(len(self.star_systems))
        for filt in self.filt_names:
            self.star_systems[filt] += rand_red * delta_red_filt[filt]

            # If ifmf is specified, then return photometry of compact objects to 0
            if self.ifmf != None:
                compact = np.where(self.star_systems['Rem_ID'] > 0)
                self.star_systems[filt][compact] = 0

        # Perturb the companions by the same amount.
        if self.imf.make_multiples:
            rand_red_comp = np.repeat(rand_red, self.star_systems['N_companions'])
            assert len(rand_red_comp) == len(self.companions)
            for filt in self.filt_names:
                self.companions[filt] += rand_red_comp * delta_red_filt[filt]

        # Finally, we'll add a column to star_systems with the overall AKs for each star
        diff_AKs = deltaAKs * rand_red
        final_AKs = AKs + diff_AKs
        col = Column(final_AKs, name='AKs_f')
        self.star_systems.add_column(col)
        #t2 = time.time()
        #print 'Diff redden: {0}'.format(t2 - t1)
        return

class ResolvedClusterDiffRedden2(ResolvedCluster):
    """
    Same as the other differentially reddened cluster, but allowing
    for asymmetric dAKs distribution
    """
    def __init__(self, iso, imf, cluster_mass, deltaAKs_blue, deltaAKs_red,
                 red_law=default_red_law, filters=None, verbose=False):

        ResolvedCluster.__init__(self, iso, imf, cluster_mass, filters=filters, verbose=verbose)

        # For a given delta_AKs (sigma of reddening distribution at Ks),
        # figure out the equivalent delta_filt values for all other filters.
        delta_red_blue_filt = {}
        delta_red_red_filt = {}
        AKs = iso.points.meta['AKS']
        red_vega_blue = vega * red_law.reddening(AKs - deltaAKs_blue).resample(vega.wave)
        red_vega = vega * red_law.reddening(AKs).resample(vega.wave)
        red_vega_red = vega * red_law.reddening(AKs + deltaAKs_red).resample(vega.wave)

        for filt in self.filt_names:
            #=======BUG HERE=====#
            # This is only because I changed column header conventions in popstar after
            # generating the IMF grid. I will need to fix this for new isochrones
            #if filt == 'mag_127m':
            #    filt_info = get_filter_info(iso.filters['hst_F127M'])
            #elif filt == 'mag_153m':
            #    filt_info = get_filter_info(iso.filters['hst_F153M'])
            #elif filt == 'mag_139m':
            #    filt_info = get_filter_info(iso.filters['hst_F139M'])
            #====================#
            filt_info = get_filter_info(iso.filters[filt.replace('mag_', '')])
            
            mag_blue = mag_in_filter(red_vega_blue, filt_info)
            mag = mag_in_filter(red_vega, filt_info)
            mag_red = mag_in_filter(red_vega_red, filt_info)
            delta_red_blue_filt[filt] = mag - mag_blue
            delta_red_red_filt[filt] = mag_red - mag

        #----Perturb all of star systems' photometry by a random amount----#
        rand_red = np.random.randn(len(self.star_systems))
        # Identify positive and negative rand_red values. This
        # will determine if the star gets a blue or red side
        # dAks, respectively.
        blue = np.where(rand_red >= 0)
        red = np.where(rand_red < 0)
        # For each star, generate a random number from 0-1. This will
        # be used to set the dAKs value through the inverse CDF
        samp = np.random.random_sample(len(self.star_systems))
        for filt in self.filt_names:
            blue_dist = scipy.stats.halfnorm(loc=0, scale=delta_red_blue_filt[filt])
            red_dist = scipy.stats.halfnorm(loc=0, scale=delta_red_red_filt[filt])

            blue_ppf = blue_dist.ppf(samp[blue])
            red_ppf = red_dist.ppf(samp[red])

            # Add appropriate dAKs value to star mag
            self.star_systems[filt][blue] -= blue_ppf
            self.star_systems[filt][red] += red_ppf
            
        # Perturb the companions by the same amount.
        #if self.imf.make_multiples:
        #    rand_red_comp = np.repeat(rand_red, self.star_systems['N_companions'])
        #    assert len(rand_red_comp) == len(self.companions)
        #    for filt in self.filt_names:
        #        self.companions[filt] += rand_red_comp * delta_red_filt[filt]
            
        return

class UnresolvedCluster(Cluster):
    def __init__(self, iso, imf, cluster_mass,
                 wave_range=[5000, 50000], verbose=False):
        """
        iso : Isochrone
        """
        # Doesn't do much.
        Cluster.__init__(self, iso, imf, cluster_mass, verbose=verbose)
        
        # Sample a power-law IMF randomly
        self.mass, isMulti, compMass, sysMass = imf.generate_cluster(cluster_mass)
        
        temp = np.zeros(len(self.mass), dtype=float)
        self.mass_all = np.zeros(len(self.mass), dtype=float)
        self.spec_list = [None] * len(self.mass)
        # placeholder array to make spectrum summing more efficient
        spec_list_np = np.zeros(shape=(len(iso.spec_list[0].flux),len(self.mass)), dtype=float)
        self.spec_list_trim = [None] * len(self.mass)
        # same as spec_list_np, but for the wavelength-trimmed spectra
        trimtmp = spectrum.trimSpectrum(iso.spec_list[0],wave_range[0],wave_range[1])
        trimx = len(trimtmp._fluxtable)
        spec_list_trim_np = np.zeros(shape=(trimx,len(self.mass)), dtype=float)

        t1 = time.time()
        for ii in range(len(self.mass)):
            # Find the closest model mass (returns None, if nothing with dm = 0.1
            mdx = match_model_mass(iso.points['mass'], self.mass[ii])
            if mdx == None:
                continue

            # getting the temp, mass, spectrum of the matched star
            temp[ii] = iso.points['Teff'][mdx]
            self.mass_all[ii] = iso.points['mass'][mdx]
            tmpspec = iso.spec_list[mdx]

            # resampling the matched spectrum to a common wavelength grid
            tmpspec = spectrum.CompositeSourceSpectrum.tabulate(tmpspec)
            tmpspecresamp = spectrum.TabularSourceSpectrum.resample(tmpspec,iso.spec_list[0].wave)
            self.spec_list[ii] = tmpspecresamp
            spec_list_np[:,ii]=np.asarray(tmpspecresamp._fluxtable)

            # and trimming to the requested wavelength range
            tmpspectrim = spectrum.trimSpectrum(tmpspecresamp,wave_range[0],wave_range[1])
            self.spec_list_trim[ii] = tmpspectrim
            spec_list_trim_np[:,ii] = np.asarray(tmpspectrim._fluxtable)
            

        t2 = time.time()
        print( 'Mass matching took {0:f} s.'.format(t2-t1))

        # Get rid of the bad ones
        idx = np.where(temp != 0)[0]
        cdx = np.where(temp == 0)[0]

        self.mass_all = self.mass_all[idx]
        self.spec_list = [self.spec_list[iidx] for iidx in idx]
        spec_list_np = spec_list_np[:,idx]
        self.spec_list_trim = [self.spec_list_trim[iidx] for iidx in idx]
        spec_list_trim_np = spec_list_trim_np[:,idx]

        self.spec_tot_full = np.sum(spec_list_np,1)

        t3 = time.time()
        print( 'Spec summing took {0:f}s'.format(t3-t2))

        self.spec_trim = np.sum(spec_list_trim_np,1)
        
        t4 = time.time()
        print( 'Spec trimming took {0:f}s'.format(t4-t3))

        self.mass_tot = np.sum(sysMass[idx])
        print( 'Total cluster mass is {0:f} M_sun'.format(self.mass_tot))

        return
        
class Isochrone(object):
    def __init__(self, logAge, AKs, distance,
                 evo_model=default_evo_model, atm_func=default_atm_func,
                 red_law=default_red_law, mass_sampling=1,
                 wave_range=[5000, 42500], min_mass=None, max_mass=None):
        """
        Parameters
        ----------
        logAge : float
            The log of the age of the isochrone.
        AKs : float
            The extinction in units if A_Ks (mag).
        distance : float
            The distance in pc.
        evModel : model cl
        mass_sampling - Sample the raw isochrone every ## steps. The default
                       is mass_sampling = 10, which takes every 10th point.
                       The isochrones are already very finely sampled. Must be
                       an integer value.
        wave_range : list
            length=2 list with the wavelength min/max of the final spectra.
            Units are Angstroms. 
        min_mass: float or None
            If float, defines the minimum mass in the isochrone.
            Units: solar masses
        max_mass: float or None
            If float, defines the maxmimum mass in the isochrone.
            Units: solar masses
        """

        t1 = time.time()
        
        c = constants

        # Get solar metallicity models for a population at a specific age.
        # Takes about 0.1 seconds.
        evol = evo_model.isochrone(age=10**logAge)  # solar metallicity

        # Eliminate cases where log g is less than 0
        idx = np.where(evol['logg'] > 0)
        evol = evol[idx]

        # Trim to desired mass range
        if min_mass != None:
            idx = np.where(evol['mass'] >= min_mass)
            evol = evol[idx]
        if max_mass != None:
            idx = np.where(evol['mass'] <= max_mass)
            evol = evol[idx] 

        # Trim down the table by selecting every Nth point where
        # N = mass sampling factor.
        evol = evol[::mass_sampling]

        # Determine which stars are WR stars.
        keys = evol.keys()
        if 'logT_WR' in keys:
            evol['isWR'] = evol['logT'] != evol['logT_WR']
            isWR_all = evol['isWR']
        else:
            isWR_all = ['None'] * len(evol)

        # Give luminosity, temperature, mass, radius units (astropy units).
        L_all = 10**evol['logL'] * c.L_sun # luminsoity in erg/s
        T_all = 10**evol['logT'] * units.K
        R_all = np.sqrt(L_all / (4.0 * math.pi * c.sigma_sb * T_all**4))
        mass_all = evol['mass'] * units.Msun # masses in solar masses
        logg_all = evol['logg']

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

        t2 = time.time()
        print( 'Isochrone generation took {0:f} s.'.format(t2-t1))
        
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
                 min_mass=None, max_mass=None, rebin=True, recomp=False, 
                 filters={'wfc3,ir,f127m', 'wfc3,ir,f139m',
                          'wfc3,ir,f153m', 'acs,wfc1,f814w',
                          'wfc3,ir,f125w', 'wfc3,ir,f160w',
                          'nirc2,J', 'nirc2,H', 'nirc2,Kp',}):

        """
        Make an isochrone with photometry in various filters.

        Description
        -----------
        Make an isochrone with photometry in various filters. Load from file
        or save to file if possible.

        Parameters
        ---------- 
        rebin: boolean (default=True)
            If true, rebins the filter functions such that they have no more than 1500 pts
            over the non-zero throughput region


        Returns
        -------
                 
        """
        
        # Make and input/output file name for the stored isochrone photometry.
        save_file_fmt = '{0}/iso_{1:.2f}_{2:4.2f}_{3:4s}.fits'
        self.save_file = save_file_fmt.format(iso_dir, logAge, AKs, str(distance).zfill(5))

        # Expected filters
        self.filters = filters

        if (not os.path.exists(self.save_file)) | (recomp==True):
            Isochrone.__init__(self, logAge, AKs, distance,
                               evo_model=evo_model, atm_func=atm_func,
                               red_law=red_law, mass_sampling=mass_sampling,
                               min_mass=min_mass, max_mass=max_mass)
            self.verbose = True
            
            # Make photometry
            self.make_photometry(rebin=rebin, vega=vega)
        else:
            self.points = Table.read(self.save_file)
            # Add some error checking.

        return

    def make_photometry(self, rebin=True, vega=vega):
        """ 
        Make synthetic photometry for the specified filters. This function
        udpates the self.points table to include new columns with the
        photometry.
        
        """
        startTime = time.time()

        meta = self.points.meta

        print( 'Making photometry for isochrone: log(t) = %.2f  AKs = %.2f  dist = %d' % \
            (meta['LOGAGE'], meta['AKS'], meta['DISTANCE']))
        print( '     Starting at: ', datetime.datetime.now(), '  Usually takes ~5 minutes')

        npoints = len(self.points)
        verbose_fmt = 'M = {0:7.3f} Msun  T = {1:5.0f} K  m_{2:s} = {3:4.2f}'

        # Loop through the filters, get filter info, make photometry for
        # all stars in this filter.
        for ii in self.filters:
            prt_fmt = 'Starting filter: {0:s}   Elapsed time: {1:.2f} seconds'
            print( prt_fmt.format(ii, time.time() - startTime))

            filt = get_filter_info(ii, rebin=rebin, vega=vega)
            filt_name = get_filter_col_name(ii)

            # Make the column to hold magnitudes in this filter. Add to points table.
            col_name = 'm_' + filt_name
            mag_col = Column(np.zeros(npoints, dtype=float), name=col_name)
            self.points.add_column(mag_col)
            
            # Loop through each star in the isochrone and do the filter integration
            print('Starting synthetic photometry')
            for ss in range(npoints):
                star = self.spec_list[ss]  # These are already extincted, observed spectra.
                star_mag = mag_in_filter(star, filt)
                
                self.points[col_name][ss] = star_mag
        
                if (self.verbose and (ss % 100) == 0):
                    print( verbose_fmt.format(self.points['mass'][ss], self.points['Teff'][ss],
                                             filt_name, star_mag))

        endTime = time.time()
        print( '      Time taken: {0:.2f} seconds'.format(endTime - startTime))

        if self.save_file != None:
            self.points.write(self.save_file, overwrite=True)

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

#===================================================#
# Iso table: same as IsochronePhot object, but doesn't do reddening application
# or photometry automatically. These are separate functions on the object.
#===================================================#
class iso_table(object):
    def __init__(self, logAge, distance, evo_model=default_evo_model,
                 atm_func=default_atm_func, mass_sampling=1,
                 min_mass=None, max_mass=None, wave_range=[5000, 30000],
                 rebin=True):
        """
        Generate an isochrone table containing star mass, temp, radius,
        luminosity, and logg, as well as a table of spectra for those
        stars. Also produce set of corresponding spectra which are
        flux calibrated by distance but not reddened in any way.

        Functions on this object:
        apply_reddening (Apply reddening with defined redlaw and AKs)
        make_photometry (make synthetic photometry for spectra)
                 
        Parameters
        ----------
        logAge : float
            The log of the age of the isochrone.
        distance : float
            The distance in pc.
        evo_model : PopStar evolution object
            Stellar evolution models used
        atm_func: PopStar atmosphere object
            Atmospheric models used
        mass_sampling - Sample the raw isochrone every ## steps. The default
                       is mass_sampling = 10, which takes every 10th point.
                       The isochrones are already very finely sampled. Must be
                       an integer value.
        min_mass: float or None
            If float, defines the minimum mass in the iso_table.
            Units: solar masses
        max_mass: float or None
            If float, defines the maxmimum mass in the iso_table.
            Units: solar masses
        wave_range : list
            length=2 list with the wavelength min/max of the final spectra.
            Units are Angstroms.
        dir_to_VISTA: string (default = './') or None
            Path to files which define the VISTA bandpasses. If None, will not
            have access to VISTA filters
        dir_to_DEC: string (default = './') or None
            Path to files which define the DECam bandpasses. If None, will not
            have access to DECam filter
        dir_to_PS1: string (default = './') or None
            Path to files which define the PS1 bandpasses. If None, will not
            have access to PS1 filters
        rebin: boolean
            If true, rebin the VISTA filter functions to match the synthetic
            spectrum. This is very useful to save computation time down the
            road.
        """
        t1 = time.time()        
        c = constants

        # Get solar metallicity models for a population at a specific age.
        # Takes about 0.1 seconds.
        evol = evo_model.isochrone(age=10**logAge)  # solar metallicity 
        
        # Eliminate cases where log g is less than 0
        idx = np.where(evol['logg'] > 0)
        evol = evol[idx]

        # Trim to desired mass range
        if min_mass != None:
            idx = np.where(evol['mass'] >= min_mass)
            evol = evol[idx]
        if max_mass != None:
            idx = np.where(evol['mass'] <= max_mass)
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
            
            # Save the final spectrum to our spec_list for later use.            
            self.spec_list.append(star)

        # Append all the meta data to the summary table.
        
        tab.meta['ATMFUNC'] = atm_func.__name__
        tab.meta['EVOMODEL'] = type(evo_model).__name__
        tab.meta['LOGAGE'] = logAge
        tab.meta['DISTANCE'] = distance
        tab.meta['WAVEMIN'] = wave_range[0]
        tab.meta['WAVEMAX'] = wave_range[1]

        self.points = tab
    
        t2 = time.time()
        print('Isochrone generation took {0:f} s.'.format(t2-t1))
        
        return

    def apply_reddening(self, AKs, extinction_law, dAKs=0, dist='uniform', dAKs_max=None):
        """
        Apply extinction to the spectra in iso_table, using the defined
        extinction law

        Parameters:
        ----------
        AKs: float
            Total extinction in AKs
            
        extinction_law: popstar extinction object
            Extinction law to be used on the spectra

        dAks: float (default = 0)
            Differential extinction to apply to star, if desired.
            Will draw reddening from Aks +/- dAks

        dAKs_max: float or None
            If not none, defines the maximum |dAKs| a star can
            have in gaussian distribution case 

        dist: string, 'uniform' or 'gaussian'
            Distribution to draw differential reddening from. If uniform,
            dAKs will cut off at Aks +/- dAKs. Otherwise, will draw
            from Gaussian of width AKs +/- dAks
            
        """
        self.AKs = np.ones(len(self.spec_list))
        # Apply reddening to each object in the spec list
        for i in range(len(self.spec_list)):
            star = self.spec_list[i]

            # Calculate reddening at extinction value using defined
            # extinction law
            if dAKs != 0:
                if dist == 'gaussian':
                    AKs_act = np.random.normal(loc=AKs, scale=dAKs)
                    # Apply dAKs_max if desired. Redo if diff > dAKs_max
                    if dAKs_max != None:
                        diff = abs(AKs_act - AKs)
                        while diff > dAKs_max:
                            print('While loop active')
                            AKs_act = np.random.normal(loc=AKs, scale=dAKs)
                            diff = abs(AKs_act - AKs)
                elif dist == 'uniform':
                    low = AKs - dAKs
                    high = AKs + dAKs
                    AKs_act = np.random.uniform(low=low, high=high)
                else:
                    print('dist {0} undefined'.format(dist))
                    return
            else:
                AKs_act = AKs

            red = extinction_law.reddening(AKs_act).resample(star.wave) 
            star *= red

            # Update the spectrum in spec list
            self.spec_list[i] = star
            self.AKs[i] = AKs_act

        # Update the table to reflect the AKs used
        self.points.meta['AKS'] = AKs

        return

    def make_photometry(self, filters, rebin=True):
        """ 
        Make synthetic photometry for the specified filters. This function
        udpates the self.points table to include new columns with the
        photometry.

        Parameters
        ----------
        filters : dictionary
            A dictionary containing the filter name (for the output columns)
            and the filter specification string that can be processed by pysynphot.
                   
        rebin: boolean
            True to rebin filter function (only used if non-zero transmission points are 
            larger than 1500 points)
 
        """
        npoints = len(self.points)

        # Loop through the filters, get filter info, make photometry for
        # all stars in this filter.
        ts = time.time()
        for filt_name, filt_str in filters.items():
            # Define filter info
            prt_fmt = 'Starting filter: {0:s}   Elapsed time: {1:.2f} seconds'
            print( prt_fmt.format(filt_name, time.time() - ts))
            filt = get_filter_info(filt_str, rebin=rebin, vega=vega)

            # Make the column to hold magnitudes in this filter. Add to points table.
            col_name = 'mag_' + filt_name
            mag_col = Column(np.zeros(npoints, dtype=float), name=col_name)
            self.points.add_column(mag_col)
            
            # Loop through each star in the isochrone and do the filter integration
            for ss in range(npoints):
                star = self.spec_list[ss]  # These are already extincted, observed spectra.
                star_mag = mag_in_filter(star, filt)
                
                self.points[col_name][ss] = star_mag
        

        endTime = time.time()
        print( '      Time taken: {0:.2f} seconds'.format(endTime - ts))

        return

def get_filter_info(name, vega=vega, rebin=True):
    """ 
    Define filter functions, setting ZP according to
    Vega spectrum. Input name is the popstar
    obs_string
    """
    tmp = name.split(',')
    filterName = tmp[-1]
        
    if name.startswith('nirc2'):
        filt = filters.get_nirc2_filt(filterName)

    elif name.startswith('2mass'):
        filt = filters.get_2mass_filt(filterName)
        
    elif name.startswith('vista'):
        filt = filters.get_vista_filt(filterName)

    elif name.startswith('decam'):
        filt = filters.get_decam_filt(filterName)

    elif name.startswith('ps1'):
        filt = filters.get_PS1_filt(filterName)

    elif name.startswith('jwst'):
        filt = filters.get_jwst_filt(filterName)

    elif name.startswith('jg'):
        filt = filters.get_Johnson_Glass_filt(filterName)
        
    elif name.startswith('nirc1'):
        filt = filters.get_nirc1_filt(filterName)
        
    elif name.startswith('ctio_osiris'):
        filt = filters.get_ctio_osiris_filt(filterName)

    else:
        filt = ObsBandpass(name)
        
        # Convert to ArraySpectralElement for resampling.
        filt = spectrum.ArraySpectralElement(filt.wave, filt.throughput,
                                             waveunits=filt.waveunits,
                                             name=filt.name)
           
    # If rebin=True, limit filter function to <=1500 wavelength points
    # over the non-zero values
    idx = np.where(filt.throughput > 0.001)[0]
    if rebin:
        if len(filt.wave[idx]) > 1500:
            new_wave = np.linspace(filt.wave[idx[0]], filt.wave[idx[-1]], 1500, dtype=float)
            filt = filt.resample(new_wave)

    # Check that vega spectrum covers the wavelength range of the filter.
    # Otherwise, throw an error
    idx = np.where(filt.throughput > 0.001)[0]
    if (min(filt.wave[idx]) < min(vega.wave)) | (max(filt.wave[idx]) > max(vega.wave)):
        raise ValueError('Vega spectrum doesnt cover filter wavelength range!')  

    vega_obs = obs.Observation(vega, filt, binset=filt.wave, force='taper')
    #vega_flux = vega_obs.binflux.sum()
    diff = np.diff(vega_obs.binwave)
    diff = np.append(diff, diff[-1])
    vega_flux = np.sum(vega_obs.binflux * diff)
    
    vega_mag = 0.03

    filt.flux0 = vega_flux
    filt.mag0 = vega_mag

    return filt

def get_filter_col_name(obs_str):
    """
    Get standard column name for synthetic photometry based on 
    the input string. The input string is expected to be an
    appropriate popstar obs_string
    """
    # How we deal with obs_string is slightly different depending
    # if it is an hst filter (and thus pysynphot syntax) or our
    # own defined filters
    tmp = obs_str.split(',')

    if len(tmp) == 3:
        filt_name = 'hst_{0}'.format(tmp[-1])
    else:
        filt_name = '{0}_{1}'.format(tmp[0], tmp[1])
        
    return filt_name

def get_obs_str(col):
    """
    Helper function to get the associated popstar obs_str given
    a column name
    """
    # Remove the trailing m_
    name = col[2:]
    
    # Define dictionary for filters
    filt_list = {'hst_f127m': 'wfc3,ir,f127m', 'hst_f139m': 'wfc3,ir,f139m', 'hst_f153m': 'wfc3,ir,f153m',
                 'hst_f814w': 'acs,wfc1,f814w', 'hst_f125w': 'wfc3,ir,f125w', 'hst_f160w': 'wfc3,ir,f160w',
                 'decam_y': 'decam,y', 'decam_i': 'decam,i', 'decam_z': 'decam,z',
                 'decam_u':'decam,u', 'decam_g':'decam,g', 'decam_r':'decam,r',
                 'vista_Y':'vista,Y', 'vista_Z':'vista,Z', 'vista_J': 'vista,J',
                 'vista_H': 'vista,H', 'vista_Ks': 'vista,Ks',
                 'ps1_z':'ps1,z', 'ps1_g':'ps1,g', 'ps1_r': 'ps1,r',
                 'ps1_i': 'ps1,i', 'ps1_y':'ps1,y',
                 'jwst_F090W': 'jwst,F090W', 'jwst_F164N': 'jwst,F164N', 'jwst_F212N': 'jwst,F212N',
                 'jwst_F323N':'jwst,F323N', 'jwst_F466N': 'jwst,F466N',
                 'nirc2_J': 'nirc2,J', 'nirc2_H': 'nirc2,H', 'nirc2_Kp': 'nirc2,Kp', 'nirc2_K': 'nirc2,K',
                 'nirc2_Lp': 'nirc2,Lp', 'nirc2_Ms': 'nirc2,Ms', 'nirc2_Hcont': 'nirc2,Hcont',
                 'nirc2_FeII': 'nirc2,FeII', 'nirc2_Brgamma': 'nirc2,Brgamma',
                 'jg_J': 'jg,J', 'jg_H': 'jg,H', 'jg_K': 'jg,K',
                 'nirc1_K':'nirc1,K', 'ctio_osiris_K': 'ctio_osirirs,K'}

    obs_str = filt_list[name]
        
    return obs_str

def rebin_spec(wave, specin, wavnew):
    """
    Helper function to rebin spectra, from Jessica's post
    on Astrobetter
    """
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs_f = obs.Observation(spec, filt, binset=wavnew, force='taper')
 
    return obs_f.binflux

def make_isochrone_grid(age_arr, AKs_arr, dist_arr, evo_model=default_evo_model,
                        atm_func=default_atm_func, redlaw = default_red_law,
                        iso_dir = './', mass_sampling=1,
                        filters=['wfc3,ir,f127m',
                                 'wfc3,ir,f139m',
                                 'wfc3,ir,f153m']):
    """
    Wrapper routine to generate a grid of isochrones of different ages,
    extinctions, and distances. 

    Parameters:
    ----------
    age_arr: array
        Array of ages to loop over, in log years

    Aks_arr: array
        Array of Aks values to loop over, in magnitudes

    dist_arr: array
        Array of distances to loop over (pc)
 
    evo_models: popstar evolution object
        Which evolution models to use

    atm_models: popstar atmospheres object
        Which atmosphere models to use

    redlaw: popstar reddening object
        Which reddening law to use

    iso_dir: str
        Directory to write the isochrones to

    mass_sampling: int
        Mass sampling of isochrone, relative to original mass sampling

    filters: dictionary
        Which filters to do the synthetic photometry on    
    """
    print( '**************************************')
    print( 'Start generating isochrones')
    print( 'Evolutionary Models adopted: {0}'.format(evo_model))
    print( 'Atmospheric Models adopted: {0}'.format(atm_func))
    print( 'Reddening Law adopted: {0}'.format(redlaw))
    print( 'Isochrone Mass sampling: {0}'.format(mass_sampling))
    print( '**************************************')

    num_models = len(age_arr) * len(AKs_arr) * len(dist_arr)
    iteration = 0
    #Loop structure: loop 1 = age, loop 2 = Aks, loop 3 = distance
    for i in range(len(age_arr)):
        for j in range(len(AKs_arr)):
            for k in range(len(dist_arr)):
                    iso = IsochronePhot(age_arr[i], AKs_arr[j], dist_arr[k],
                                        evo_model=evo_model, atm_func=atm_func,
                                        red_law=redlaw, iso_dir=iso_dir,
                                        mass_sampling=mass_sampling,
                                        filters=filters)
                    iteration += 1
                    print( 'Done ' + str(iteration) + ' of ' + str(num_models))

    # Also, save a README file in iso directory documenting the params used
    _out = open(iso_dir+'README.txt', 'w')
    _out.write('Popstar parameters used to generate isochrone grid:\n')
    _out.write('Evolutionary Models: {0}\n'.format(evo_model))
    _out.write('Atmospheric Models: {0}\n'.format(atm_func))
    _out.write('Reddening Law: {0}\n'.format(redlaw))
    _out.write('Isochrone Mass: {0}'.format(mass_sampling))
    _out.close()
    return

# Little helper utility to get the magnitude of an object through a filter.
def mag_in_filter(star, filt):
    """
    Assumes that extinction is already resampled to same wavelengths
    as filter, and has been applied.
    """
    star_in_filter = obs.Observation(star, filt, binset=filt.wave, force='taper')
    #star_flux = star_in_filter.binflux.sum()
    diff = np.diff(star_in_filter.binwave)
    diff = np.append(diff, diff[-1])
    star_flux = np.sum(star_in_filter.binflux * diff)
    
    star_mag = -2.5 * math.log10(star_flux / filt.flux0) + filt.mag0
    return star_mag

def match_model_mass(isoMasses,theMass):
    dm = np.abs(isoMasses - theMass)
    mdx = dm.argmin()

    # Model mass has to be within 2% of the desired mass
    if (dm[mdx] / theMass) > 0.1:
        return None
    else:
        return mdx

def match_model_masses(isoMasses, starMasses):
    kdt = KDTree( isoMasses.reshape((len(isoMasses), 1)) )
    q_results = kdt.query(starMasses.reshape((len(starMasses), 1)), k=1)
    indices = q_results[1]

    dm_frac = np.abs(starMasses - isoMasses[indices]) / starMasses

    idx = np.where(dm_frac > 0.1)[0]
    indices[idx] = -1
    
    return indices

    
def get_evo_model_by_string(evo_model_string):
    return getattr(evolution, evo_model_string)


