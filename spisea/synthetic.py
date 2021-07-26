import numpy as np
import pylab as plt
from spisea import reddening
from spisea import evolution
from spisea import atmospheres as atm
from spisea import filters
from spisea.imf import imf, multiplicity
from scipy import interpolate
from scipy import stats
from scipy.special import erf
from pysynphot import spectrum
from pysynphot import ObsBandpass
from pysynphot import observation as obs
import pysynphot
from astropy import units
from astropy.table import Table, Column, MaskedColumn, hstack, vstack
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
from scipy.spatial import cKDTree as KDTree
import inspect
import astropy.modeling
from astropy import constants
import pdb
import copy
from scipy.interpolate import LinearNDInterpolator

cs = constants
un = units
c = constants
default_evo_model = evolution.MISTv1()
default_red_law = reddening.RedLawNishiyama09()
default_atm_func = atm.get_merged_atmosphere
default_wd_atm_func = atm.get_wd_atmosphere

    
def Vega():
    
    
    # Use Vega as our zeropoint... assume V=0.03 mag and all colors = 0.0
    # These parameters are defined in Girardi+02
    vega = atm.get_kurucz_atmosphere(temperature=9550, 
                                     gravity=3.95,
                                     metallicity=-0.5)

    # Following the K93 README, set wavelength range to 0.1 - 10 microns.
    # This defines the maximum allowed wavelength range in SPISEA
    vega = spectrum.trimSpectrum(vega, 995, 100200)

    # This is (R/d)**2 as reported by Girardi et al. 2002, page 198, col 1.
    # and is used to convert to flux observed at Earth.
    vega *= 6.247e-17 
    
    return vega

vega = Vega()

class BypasserContainer():
    description = "This is so vectorization doesn't think that I'm an array"

    

class Cluster(object):
    
    
    """
    Base class to create a cluster with user-specified isochrone,
    imf, ifmr, and total mass. 

    Parameters
    -----------
    iso: isochrone object
        SPISEA isochrone object
    
    imf: imf object
        SPISEA IMF object

    cluster_mass: float
        Total initial mass of the cluster, in M_sun

    ifmr: ifmr object or None
        If ifmr object is defined, will create compact remnants
        produced by the cluster at the given isochrone age. Otherwise,
        no compact remnants are produced.

    seed: int
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output.
        Default None

    vebose: boolean
        True for verbose output.
    """
    
    
    def __init__(self, iso, imf, cluster_mass, ifmr=None, verbose=False,
                     seed=None):
        
        
        self.verbose = verbose
        self.iso = iso
        self.imf = imf
        self.ifmr = ifmr
        self.cluster_mass = cluster_mass
        self.seed = seed
        
        return
    
    
class Cluster_w_Binaries(Cluster):
    """
    Cluster sub-class that produces a *resolved* stellar cluster with
    binary evolution accounted for.
    A table is output with the synthetic photometry and intrinsic
    properties of the individual stars (or stellar systems, if
    mutliplicity is used in the IMF object).

    A second table is produced that
    contains the properties of the companion stars independent of their
    primary stars.

    Parameters
    -----------
    iso: isochrone_binary object
        SPISEA binary_isochrone object
    imf: imf object
        SPISEA IMF object
    cluster_mass: float
        Total initial mass of the cluster, in M_sun
    ifmr: ifmr object or None
        If ifmr object is defined, will create compact remnants
        produced by the cluster at the given isochrone age. Otherwise,
        no compact remnants are produced.
    seed: int
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output.
        Default None
    vebose: boolean
        True for verbose output.

    ====================
    
    """
    def __init__(self, iso, imf, cluster_mass, ifmr=None, verbose=False,
                 seed=None, tests='mass_current'):
        Cluster.__init__(self, iso, imf, cluster_mass,
                         ifmr=ifmr, verbose=verbose,
                         seed=seed)
        # Provide a user warning is random seed is set
        if seed is not None:
            print('WARNING: random seed set to %i' % seed)
        t1 = time.time()
        self.iso = iso
        self.verbose = verbose
        #####
        # Sample the IMF to build up our cluster mass.
        # NOTE: I have adjustment factor in place
        # in order to account for the amount of stars and star mass that
        # may be weeded out due to making initial masses come from the
        # isochrone
        #####

        mass, isMulti, compMass, sysMass = imf.generate_cluster(cluster_mass,seed=seed)
        self.filt_names = self.iso.filt_names_table
        # Below: INTENDED cluster mass!
        self.cluster_mass = cluster_mass
        # A little sanity check inserted to make sure that IMF generated
        # cluster mass is close to the desired cluster mass
        np_min_mass = np.min([sysMass.sum(), cluster_mass])
        #####
        single_star_systems = self.make_singles_systems_table(isMulti, sysMass)
        if self.imf.make_multiples:
            # Temporary companion mass finder:
            # Double systems is a bit of a misnomer
            # Correct term would be "multi_star_system"
            companions, double_systems = \
            self.make_primaries_and_companions(sysMass, compMass)
            self.star_systems = vstack([double_systems, single_star_systems])
            self.companions = companions
            
        else:
            self.star_systems = single_star_systems
            self.intended_primaries_mass = 0
            self.intended_companions_mass = 0
            self.intended_singles_mass = sysMass.sum()
        return

    def set_columns_of_table(self, t, N_systems, multi=False):
        """
        - Inputs -
        t:Astropy Table with entries for initial masses
        of the stars.
        N_systems: Length of AstroPy Table
        multi: whether t is the table for multi-star systems
        - Output: None
        ====================
        adds columns for basic parameters to table t with their
        descriptions (e.g. metallicity or Teff)
        """
        t.add_column(Column(np.zeros(N_systems, dtype=float),
                            name='Teff'))
        t.add_column(Column(np.empty(N_systems, dtype=float),
                            name='L'))
        t.add_column(Column(np.empty(N_systems, dtype=float),
                            name='logg'))
        t.add_column(Column(np.repeat(0.0, N_systems),
                            name='isWR'))
        t.add_column(Column(np.empty(N_systems, dtype=float),
                            name='mass_current'))
        t.add_column(Column(np.empty(N_systems, dtype=float),
                            name='phase'))
        t.add_column(Column(np.repeat(False, N_systems),
                            name='IFMR_it'))
        # To ADD: 
        # t.add_column(Column(np.repeat(False, N_systems),
        #.                    name='out_of_range'))
        t.add_column(Column(np.repeat(self.iso.metallicity,
                                      N_systems),
                            name='metallicity'))
        t.add_column(Column(np.repeat(multi, N_systems),
                                       name='isMultiple'))
        t.add_column(Column(np.repeat(0.0, N_systems),
                            name='merged'))
        # Add the filter columns to the table. They are empty so far.
        # Keep track of the filter names in : filt_names
        for filt in self.filt_names:
            t.add_column(Column(np.zeros(N_systems,
                                         dtype=float),
                                name=filt))
        if multi:
            t.add_column(Column(np.arange(N_systems),
                                name='designation'))
            
        return None

    def applying_IFMR_stars(self, stars):
        """
        ====================
        Applies the IFMR (if there's any given) to the
        a table of stars. Specifically, when stars have phase of 5 but
        have an effective temperature of 0 or some non-finite quantity,
        or have logg >= 6.9 or have non-physical temperatures
        (and are not merged with the primary)
        This is in order to make sure there are no compact stars masquerading
        as non-compact stars.

        Input: stars: Astropy table of stars (either of single stars, primary stars,
        or companions)

        Output: None
        """
        star_systems_phase_non_nan = np.nan_to_num(stars['phase'],
                                                   nan=-99)
        bad = np.where((star_systems_phase_non_nan > 5) &
                       (star_systems_phase_non_nan < 101) &
                       (star_systems_phase_non_nan != -99))
        # Print warning, if desired
        if self.verbose:
            for ii in range(len(bad[0])):
                print('WARNING: changing phase ' +
                      '{0} to 5'.format(star_systems_phase_non_nan[bad[0]
                                        [ii]]))
        stars['phase'][bad] = 5
        # Also dealing with interpolation screwing up phases less than 5
        bad2 = np.where((star_systems_phase_non_nan < 5))
        # Print warning, if desired
        if self.verbose:
            for ii in range(len(bad[0])):
                print('WARNING: changing phase ' +
                      '{0} to 5'.format(star_systems_phase_non_nan[bad[0]
                                        [ii]]))
        stars['phase'][bad2] = -99
        ####
        # The IFMR is applied on stars with phase 5 (which are therefore not
        # secondary stars that were merged with the system primary) and have
        # have nonphysical effective temperature (Teff = 0 or some finite value)
        # and/or have surface gravity of at least 6.9 cgs.
        #####
        
        cdx_rem = np.where(((stars['phase'] == 5) & ((stars['Teff'] == 0) |
                                                     (~np.isfinite(stars['Teff'])) |
                                                     (stars['logg'] >= 6.9))) |
                           (stars['IFMR_it']))[0]
        if self.ifmr:
            # Identify compact objects as those with Teff = 0.
            # Conditions the star has to be not merged and the star has to be
            # Calculate remnant mass and ID for compact
            # objects; update remnant_id and
            # remnant_mass arrays accordingly
            cond = ('metallicity_array' in
                    inspect.getfullargspec(self.ifmr.generate_death_mass).args)
            if cond:
                r_mass_tmp, r_id_tmp = \
                self.ifmr.generate_death_mass(mass_array=stars
                                              ['mass'][cdx_rem],
                                              metallicity_array=stars
                                              ['metallicity'][cdx_rem])
            else:
                r_mass_tmp, r_id_tmp = \
                self.ifmr.generate_death_mass(mass_array=stars
                                              ['mass'][cdx_rem])
            # Drop remnants where it is not relevant
            # (e.g. not a compact object or
            # outside mass range IFMR is defined for)
            good = np.where(r_id_tmp > 0)
            cdx_rem_good = cdx_rem[good]
            stars['mass_current'][cdx_rem_good] = r_mass_tmp[good]
            # when we have somee undefined or 0 temperature,
            # shouldn't we practically have no light coming
            # out of it (not even blackbody radiation)?
            stars['L'][cdx_rem] = 0.0
            stars['phase'][cdx_rem_good] = r_id_tmp[good]
            for filt in self.filt_names:
                stars[filt][cdx_rem] = \
                np.full(len(cdx_rem), np.nan)
        else:
            # If the IFMR doesn't exist, then we may need to
            # get rid of the stars that would warrant usage of the IFMR.
            stars.remove_rows(cdx_rem)
        return None
                
    def generate_2body_parameters(self, star_systemsPrime, companions):
        """
        Input:star_systemsPrime:Astropy table of the primary stars of the cluster
        companions: Astropy table of the companion stars of the cluster
        Generate the log_separation, e, i, omegas,
        (2-body problem parameters) using the given
        multiplicity.
        """

        N_comp_tot = len(companions)
        inst = isinstance(self.imf._multi_props,
                          multiplicity.MultiplicityResolvedDK)
        if inst:
            companions.add_column(Column(np.zeros(N_comp_tot, dtype=float),
                                         name='log_a'))
            companions.add_column(Column(np.zeros(N_comp_tot, dtype=float),
                                         name='e'))
            companions.add_column(Column(np.zeros(N_comp_tot, dtype=float),
                                         name='i', description='degrees'))
            companions.add_column(Column(np.zeros(N_comp_tot, dtype=float),
                                         name='Omega'))
            companions.add_column(Column(np.zeros(N_comp_tot, dtype=float),
                                         name='omega'))
            for ii in range(len(companions)):
                # If the mass of the star is less than 0.1 I will
                # then let the log_a be 0.
                # Do not want to trigger an error but I will try
                # NOT to cause and instead make
                # log_a = np.nan. That will be our indicator that
                # I could not calculate the log_a using
                # Duchene and Krauss and that we will
                # not use the log(a)
                try:
                    companions['log_a'][ii] = \
                    self.imf._multi_props.log_semimajoraxis(star_systemsPrime
                                                            ['mass']
                                                            [companions
                                                             ['system_idx']
                                                             [ii]])
                except ValueError:
                    # Indicator that we can't use the Multiplicity for log_a
                    companions['log_a'][ii] = np.nan
                    continue
            companions['e'] = \
            self.imf._multi_props.random_e(np.random.rand(N_comp_tot))
            props = self.imf._multi_props
            companions['i'], companions['Omega'], companions['omega'] = \
            props.random_keplarian_parameters(np.random.rand(N_comp_tot),
                                              np.random.rand(N_comp_tot),
                                              np.random.rand(N_comp_tot))
        else:
            # Indicator that we can't use the Multiplicity for log_a
            companions['log_a'] = np.nan
            companions['e'] = np.nan
            companions['i'], companions['Omega'], companions['omega'] = \
            np.nan, np.nan, np.nan
        return None
    
    def finding_secondary_stars(self, star_systemsPrime, companions):
        """
        ====================
        Given the table so far of the primary stars and of the companions
        we try to find the secondary of the primary star. Secondary star
        is the companion star with a) greatest gravitational pull on the primary
        (if isochrone provides primary) or b) greatest mass.

        Input: star_systemsPrime: table of primary stars (with companions)
        companions: table of companions (currently including corresponding index
        
        Output: compMass_IDXs: a hashmap linking the star system to the next companion
        star in compMass array (generated by IMF) that needs to be evaluated in
        filling_in_primaries_and_companions
        max_log_gs: hashmap linking star system to a 2-element list such that the first
        element of the list is the GM_secondary star/(separation)^2 and Mass of the secondary star
        
        
        """
        compMass_IDXs = {}
        max_log_gs = {}
        
        for x in range(len(star_systemsPrime)):
            # subtbl stands for sub-table. (of companions whose primary star
            # is at index x.
            subtbl = companions[np.where(companions['system_idx'] == x)[0]]
            subtbl2 = copy.deepcopy(subtbl)
            max_log_gs[x] = [0, 0]
            # handling nan'ed separations. We generally don't want to
            # consider for minimum distance and maximum gravitational
            # influence to the primary
            subtbl2['log_a'][np.where(~np.isfinite(subtbl['log_a']))] = np.inf
            max_log_gs[x][0] = np.max(subtbl2['mass']/(10 ** subtbl2['log_a']) ** 2)
            max_log_gs[x][1] = np.max(subtbl2['mass'])
              
            # Indicate where in the compMass array for star system
            # number x we are in
            # when we start looking at the companion
            compMass_IDXs[x] = 0
            # This is where the matching of the primary and the secondary begin
            # I will be tracking the number of bad star systems throguh the
            # variable rejected_system and rejected_companion
        return compMass_IDXs, max_log_gs

    def filling_in_primaries_and_companions(self, star_systemsPrime,
                                            companions, compMass_IDXs,
                                            min_log_gs, compMass,
                                            error_check='mass_current'):
        """
        ====================
        Utilizes linear interpolation to try to find the luminosity, current mass,
        photometry, phase, surface gravity, and temperature of the cluster.

        Input: star_systems: AstroPy Table containing information so far
        about IMF generated star systems and their primary stars 
        companions: Table containing information so far about IMF
        generated companions.

        Output: rejected_system[:, 0], rejected_system[:, 1]: NumPy table
        The first two being the imf generated masses of primary stars
        and secondary stars such that when mass_current
        is interpolated over  the primary star's initial mass and the secondary star's
        initial mass, the primary star's or secondary star's error_check parameter is NOT
        interpolated correctly (i.e. spits out an NaN). This does not include
        stars that interpolation predicts to be the 
        good_system: NumPy table where each row contains the primary-secondary star initial
        mass pair (IMF generated) of star systems not included in rejected_system
        """

        # Case where multiplicity DK is used.
        li = list(zip(self.iso.primaries['mass'], self.iso.secondaries['mass'],
                      self.iso.secondaries['log_a']))
        li2 = list(zip(self.iso.primaries['mass'],
                       self.iso.secondaries['mass']))
        interp_keys = ['Teff', 'L', 'logg', 'isWR', 'mass_current', 'phase'] + self.filt_names
        iso_interpsP = {}
        iso_interpsS = {}
        interp_creation_time = []
        inst = isinstance(self.imf._multi_props,
                          multiplicity.MultiplicityResolvedDK)
        if not inst:
            for ikey in interp_keys:
                print("Creating Interpolator for: " + ikey)
                iso_interpsP[ikey] = LinearNDInterpolator(li2, self.iso.primaries[ikey])
                iso_interpsS[ikey] = LinearNDInterpolator(li2, self.iso.secondaries[ikey])
            iso_interpsS['log_a'] = LinearNDInterpolator(li2, self.iso.secondaries['log_a'])
            iso_interpsS['merged'] = LinearNDInterpolator(li2, self.iso.secondaries['merged'])
        else:
            for ikey in interp_keys:
                print("Creating Interpolator for: " + ikey)
                iso_interpsP[ikey] = LinearNDInterpolator(li, self.iso.primaries[ikey])
                iso_interpsS[ikey] = LinearNDInterpolator(li, self.iso.secondaries[ikey])
            iso_interpsS['merged'] = LinearNDInterpolator(li, self.iso.secondaries['merged'])
        # Dealing with single stars and higher order companions
        iso_interps3 = {}
        interp_creation_time3 = []
        for ikey in interp_keys:
            iso_interps3[ikey] = interpolate.interp1d(self.iso.singles['mass'],
                                                      self.iso.singles[ikey],
                                                      kind='linear', bounds_error=False,
                                                      fill_value=np.nan)
        # Now time to introduce it to the IMF generated stars
        rejected_system = []
        good_system = []
        if inst:
            for sys in range(len(star_systemsPrime)):
                # segundaria: index of the secondary
                comps_of_sys= np.where(companions['system_idx'] == star_systemsPrime[sys]['designation'])[0]
                segundaria = np.where((companions['system_idx'] == star_systemsPrime[sys]['designation']) &
                                      np.isclose(companions['mass'] /
                                                 (10 ** (2 * companions['log_a'])),
                                                 min_log_gs[sys][0], rtol=1e-2))[0][0]
                for k in interp_keys:
                    # Run the interpolator on all companions. Then, rerun the interpolator on the
                    # secondary star. Although it may seem redundant. Note that there is ony one
                    # secondary star and that this prevents me from having to do some really
                    # crafty and nasty logic trick.
                    companions[k][comps_of_sys] = np.round(iso_interps3[k](companions['mass'][comps_of_sys]))
                    put_in_table = iso_interpsP[k](star_systemsPrime[sys]['mass'],
                                                        companions['mass'][segundaria],
                                                        companions['log_a'][segundaria])
                    put_in_table2 = iso_interpsS[k](star_systemsPrime[sys]['mass'],
                                                        companions['mass'][segundaria],
                                                        companions['log_a'][segundaria])
                    if (k == 'phase' or k == 'isWR'):                        
                        star_systemsPrime[k][sys] = np.round(put_in_table)
                        companions[k][segundaria] = np.round(put_in_table2)
                    else:
                        star_systemsPrime[k][sys] = put_in_table
                        companions[k][segundaria] = put_in_table2
                         
                companions['merged'][segundaria] = np.round(iso_interpsS['merged']
                                                            (star_systemsPrime[sys]['mass'],
                                                             companions['mass'][segundaria],
                                                             companions['log_a'][segundaria]))
                companions['the_secondary_star?'][segundaria] = True
                if (error_check == 'Teff' or error_check == 'L'):
                    # We need to make sure that the primary did not interpolate
                    # correctly or the secondary did not interpolate correctly
                    # Consider where the star system did not interpolate correctly.
                    # Note that current parameters can be validly equal to NaN if they are of
                    # a merged secondary
                    cond_bad = (np.isnan(companions[error_check][segundaria]) and 
                                (companions['phase'][segundaria] < 101) and
                                (companions['merged'][segundaria] == 0))
                    cond_bad = cond_bad or ((np.isnan(star_systemsPrime[error_check][sys]) and
                                             (star_systemsPrime['phase'][sys] < 101)))
                else:
                    # Consider where the star system did not interpolate correctly.
                    # Note that current parameters can be validly equal to NaN if they are of
                    # a merged secondary
                    cond_bad = ((np.isnan(companions[error_check][segundaria]) and
                                 (companions['merged'][segundaria] == 0)) or
                                np.isnan(star_systemsPrime[error_check][sys]))
                if cond_bad:
                    rejected_system.append([star_systemsPrime['mass'][sys],
                                            companions['mass'][segundaria]])
                else:
                    good_system.append([star_systemsPrime['mass'][sys],
                                        companions['mass'][segundaria]])
        else:
            # Case where the multiplicity does not create the log-
            for sys in range(len(star_systemsPrime)):
                # segundaria: index of the secondary
                comps_of_sys= np.where(companions['system_idx'] == star_systemsPrime[sys]['designation'])[0]
                segundaria = np.where((companions['system_idx'] == star_systemsPrime[sys]['designation']) &
                                       np.isclose(companions['mass'],
                                                  min_log_gs[sys][1], rtol=1e-2))[0][0]
                for k in interp_keys:
                    companions[k][comps_of_sys] = iso_interps3[k](companions['mass'][comps_of_sys])
                    put_in_table = iso_interpsP[k](star_systemsPrime[sys]['mass'],
                                                    companions['mass'][segundaria])
                    put_in_table2 = iso_interpsS[k](star_systemsPrime[sys]['mass'],
                                                    companions['mass'][segundaria])
                    # we want to make sure that the phase and status of whether a star
                    # isWR are integers, not some wacky float
                    if k == 'phase' or k == 'isWR':
                        put_in_table = np.round(put_in_table)
                        put_in_table2 = np.round(put_in_table2)
                    star_systemsPrime[k][sys] = put_in_table
                    companions[k][sys] = put_in_table2
                companions['log_a'][segundaria] = iso_interpsS['log_a'](star_systemsPrime[sys]['mass'],
                                                                        companions['mass'][segundaria])
                # Check if the secondary star is merged.
                companions['merged'][segundaria] = np.round(iso_interpsS['merged']
                                                            (star_systemsPrime[sys]['mass'],
                                                             companions['mass'][segundaria]))
                # Kind of obvious, but we want to mark the secondary star as the 
                # secondary; otherwise, the star is just another companion.
                companions['the_secondary_star?'][segundaria] = True
                # Consider stars that I suspect are NOT properly inteprolated
                if (error_check == 'Teff' or error_check == 'L'):
                    # We need to make sure that the primary did not interpolate
                    # correctly or the secondary did not interpolate correctly
                    # Consider where the star system did not interpolate correctly.
                    # Note that current parameters can be validly equal to NaN if they are of
                    # a merged secondary
                    cond_bad = (np.isnan(companions[error_check][segundaria]) and 
                                (companions['phase'][segundaria] < 101) and
                                (companions['merged'][segundaria] == 0))
                    cond_bad = cond_bad or (np.isnan(star_systemsPrime[error_check][sys]) and
                                            (star_systemsPrime['phase'][sys] < 101))
                else:
                    cond_bad = ((np.isnan(companions[error_check][segundaria]) and
                                 (companions['merged'][segundaria] == 0)) or
                                np.isnan(star_systemsPrime[error_check][sys]))
                if cond_bad:
                    # If not already, note that current parameters can be validly equal to
                    # NaN if they are of a merged secondary
                    rejected_system.append([star_systemsPrime['mass'][sys],
                                            companions['mass'][segundaria]])
                else:
                    good_system.append([star_systemsPrime['mass'][sys],
                                        companions['mass'][segundaria]])
        star_systemsPrime['metallicity'] = np.ones(len(star_systemsPrime)) * self.iso.metallicity
        companions['metallicity'] = np.ones(len(companions)) * self.iso.metallicity
        rejected_system = np.array(rejected_system)
        good_system = np.array(good_system)
        return rejected_system[:, 0], rejected_system[:, 1], good_system

    def adding_up_photometry(self, star_systemsPrime, companions):
        """
        - Input - 
        star_systemsPrime: AstroPy table representing non-single
        star systems/primaries filled in partially by
        filling_in_primaries_and_companions
        companions: AstroPy table representing companion stars of
        systems represented in star_systemsPrime
        filled in partially by filling_in_primaries_and_companons
        - Output -
        None
        ====================
        Adds up initial masses of secondaries to initial mass
        of the primary to create the systemMass values.
        Also adding up the photometries (kind of hard and ugly)
        to find each star systems' aggregate photometry.
        """

        for x in range(len(star_systemsPrime)):
            # Find all companions of the star system
            sub_tbl = companions[np.where(companions['system_idx'] == x)[0]]
            sum_of_comp = sub_tbl['mass'].sum()
            # Obtain the system mass.
            star_systemsPrime['systemMass'][x] = sum_of_comp + \
            star_systemsPrime['mass'][x]
            if np.isnan(star_systemsPrime['phase'][x]):
                star_systems_phase_non_nan = -99
            else:
                star_systems_phase_non_nan = star_systemsPrime['phase'][x]
            cond = ((int(star_systems_phase_non_nan) > 5) and
            (int(star_systems_phase_non_nan) < 101) and
            (int(star_systems_phase_non_nan) != -99))
            if (cond):
                print("Changing phase of primaries")
                if (self.verbose):
                    print('WARNING: changing phase' +
                          ' {0} to 5'.format(star_systems_phase_non_nan))
                star_systemsPrime['phase'][x] = 5
            for filt in self.filt_names:
                # Magnitude of primary star (initially)
                mag_s = star_systemsPrime[filt][x]
                # Companion stars corresponding to the primary star
                # Note that I try tio make sure that 
                comps = companions[np.where((companions['system_idx'] == x) &
                                            np.isfinite(companions[filt]))[0]]
                # trying to obtain as many finite magnitudes as I can
                if (not len(comps[filt])):
                    mag_c = np.nan
                else:
                    # Finding the sum of fluxes and then taking magnitude.
                    mag_c = comps[filt]
                # Add companion flux to system flux.
                f1 = 10 ** (-1 * mag_s / 2.5)
                f2 = 10 ** (-1 * mag_c / 2.5)

                f1 = np.nan_to_num(f1)
                f2 = np.nan_to_num(f2)
                f2 = np.sum(f2)
                # Good and bad systems
                
                # If *both* objects are dark, then keep the magnitude
                # as np.nan. Otherwise, add fluxes together
                if (f1 != 0 or f2 != 0):
                    star_systemsPrime[filt][x] = -2.5 * np.log10(f1 + f2)
                else:
                    star_systemsPrime[filt][x] = np.nan
                

    def make_singles_systems_table(self, isMulti, sysMass):
        """
        Make a part of the star_systems table and get synthetic photometry
        for each single star system. 
        Input: isMulti: Whether a star system is a
        multi-star system. (output of the IMF) 
        sysMass: System masses of all stellar systems (generated by IMF)
        Output: A part of what becomes the self.star_systems table that contains all
        information (mass, luminosity, photometry) regarding single star systems.
        Data values are matched from rows of the self.isochrone
        based on the row's initial primary mass.
        """
        # We will be only looking for stars that are not multiple systems
        sysMass = sysMass[np.where(~ isMulti)[0]]
        N_systems = len(sysMass)
        star_systems = Table([sysMass, sysMass],
                             names=['mass', 'systemMass'])
        self.set_columns_of_table(star_systems, N_systems)
        # Set up for interpolation
        interp_keys = (['Teff', 'L', 'logg', 'isWR', 'mass_current', 'phase'] +
                       self.filt_names)
        iso_interps3 = {}
        for ikey in interp_keys:
            iso_interps3[ikey] = interpolate.interp1d(self.iso.singles['mass'],
                                                      self.iso.singles[ikey],
                                                      kind='linear',
                                                      bounds_error=False,
                                                      fill_value=np.nan)
        # Add columns for the Teff, L, logg, isWR,
        # mass_current, phase for the primary stars.
        star_systems['Teff'] = iso_interps3['Teff'](star_systems['mass'])
        star_systems['L'] = iso_interps3['L'](star_systems['mass'])
        star_systems['logg'] = iso_interps3['logg'](star_systems['mass'])
        star_systems['isWR'] = np.round(iso_interps3['isWR'](star_systems['mass']))
        star_systems['mass_current'] = iso_interps3['mass_current'](star_systems['mass'])
        star_systems['phase'] = np.round(iso_interps3['phase'](star_systems['mass']))
        star_systems['metallicity'] = np.ones(N_systems) * \
        self.iso.metallicity
        # Remove star systems that are outside of max mass of isochrone
        # Still inquiring about the secondary stars
        star_systems,  dummy = self._remove_bad_systems(star_systems)
        self.applying_IFMR_stars(star_systems)
        for filt in self.filt_names:
            star_systems[filt] = iso_interps3[filt](star_systems['mass'])
        # star_systems.remove_columns(['touchedP'])
        return star_systems

    def make_primaries_and_companions(self, star_systems, compMass, test='mass_current'):
        """
        - Input -
        star_systems: numpy array of all star_system
        masses (generated by IMF)
        compMass: list of lists of masses of companions of each star system

        - Output -
        Creates tables star_systemsPrime and companions,
        which contain data regarding
        star systems with multiple stars
        and the companions.

        ====================
        Makes the primary stars and companions for non-single star systems of
        the star_systems table and get synthetic
        photometry for each single star system. Creates remnants when necessary.
        Given an initial primary star mass, initial companion mass, and log
        separation (last one is included if applicable) between the star and
        the companion, we designate the star with the secondary as the star with
        least separation from the primary. If there are multiple stars with the
        same value for separation, the star with the most mass is selected out of
        the set of companions with the same minimum separation and is designated
        as the secondary.
        
        Primaries and secondaries are matched to the pairs of stars with initial
        primary mass - initial secondary mass - log_current separation from the
        isochrone that are closest to the values generated by the IMF. Given that
        the error between the isochrone star system and the imf generated system
        is small enough, the star system is included in the cluster’s star_systems
        and companions tables.
        
        Tertiary and higher order companions are included, but are matched to single
        star models. As  similar to the policy with matching the primary-secondary pairs,
        if the companion deviates too much, in this case in terms of initial mass, from
        the most similar single star from the isochrone, it does not become part of
        the cluster object.
        """

        # Obtain the indices of the systems corresponding to each companion
        # make star_systems contain only multi-star systems
        indices = [x for x in range(len(compMass)) if len(compMass[x])]
        star_systems = star_systems[indices]
        # For each star system, the total mass of companions
        compMass_sum = np.array([sum(compMass[x]) for x in indices])
        # Converting companion masses 2D list into array.
        compMass = np.array([compMass[x] for x in indices])
        # Make star_systems array only contain the masses of the primary stars
        star_systems = star_systems - compMass_sum
        # Number of multi-body star systems
        N_systems = len(star_systems)
        #####
        # MULTIPLICITY
        # Make a second table containing all the companion-star masses.
        # This table will be much longer... here are the arrays:
        # sysIndex - the index of the system this star belongs too
        # mass - the mass of this individual star.
        N_companions = np.array([len(star_masses) for star_masses in compMass])
        N_comp_tot = N_companions.sum()
        system_index = np.repeat(range(len(indices)), N_companions)
        # From now on, we will try to use the
        # new indices (i.e. what the 0-index system is
        # after we get rid of the single systems from the
        # star systems table) as much as possible.
        # Shows which star system number (index of primary star)
        # each companion star corresponds to
        companions = Table([system_index], names=['system_idx'])
        star_systemsPrime = Table([star_systems, np.zeros(N_systems,
                                                           dtype=float)],
                                  names=['mass', 'systemMass'])
        # Create a table for primary stars
        N_systems = len(indices)
        # Add columns for the Teff, L, logg, isWR,
        # mass_current, phase for the primary stars.
        self.set_columns_of_table(star_systemsPrime, N_systems, multi=True)
        
        # The following column indicates whether a star system
        # will have to be deleted due to lack of good match.
        star_systemsPrime.add_column(Column(np.repeat(False, N_systems),
                                            name='bad_system'))
        # Designation is really the index of the star,
        # but this may not be equal to
        # the 0, 1.., len -1 index iterable of the
        # cluster once we kill off the bad systems (systems whose
        # primary-secondary pair could not be found
        
        # Add columns for the Teff, L, logg, isWR mass_current,
        # phase, and filters for the companion stars.
        companions.add_column(Column(np.zeros(N_comp_tot, dtype=float),
                                     name='mass'))
        self.set_columns_of_table(companions, N_comp_tot)
        
        # Marks whether the system was unmatchable and
        # its row must be deleted.
        companions.add_column(Column(np.repeat(False, N_comp_tot),
                                     name='bad_system'))
        companions.add_column(Column(np.repeat(False, N_comp_tot),
                                     name='the_secondary_star?'))
        for ind in range(len(star_systemsPrime)):
            (companions['mass'][np.where
                                (companions['system_idx'] ==
                                 ind)]) = compMass[ind]
        
        self.generate_2body_parameters(star_systemsPrime, companions)
        # What-th companion of the system am I inspecting
        # Index of the companion we are inspecting in the compMass
        # 'Tuple' of minimum log_a, initial mass of star.
        # i.e. initial log_a, mass of secondary star of a system
        compMass_IDXs, max_log_gs = (self.finding_secondary_stars(star_systemsPrime,
                                                                  companions))
        
        self.rejected_prims, self.rejected_sec, self.good_systems = \
        self.filling_in_primaries_and_companions(star_systemsPrime,
                                                 companions, compMass_IDXs,
                                                 max_log_gs, compMass,
                                                 error_check=test)
        star_systemsPrime, companions = self._remove_bad_systems(star_systemsPrime, companions)
        
        # We don't need to worry about stars that did not get matched
        # to a close enough equivalent in the KD Tree.
        # Why? We've moved on to using the 
        
        for x in range(len(companions)):
            companions['system_idx'][x] = \
            np.where(star_systemsPrime['designation'] ==
                     companions['system_idx'][x])[0]      
        #####
        # Make Remnants for non-merged stars with nan
        # and for stars with 0 Kelvin Teff
        #####
        # Identify compact objects as those with Teff = 0 or 
        # the secondary stars that are non-merged and have a non-
        # finite temperature
        # Before applying IFMR, I mark up which stars need to be
        # should get kicked out based on the 
        self.applying_IFMR_stars(companions)
        self.applying_IFMR_stars(star_systemsPrime)
        # The compact remnants have photometric
        # fluxes of nan now. So now we can procede
        # with the photometry.
        self.adding_up_photometry(star_systemsPrime, companions)
        companions_phase_non_nan = np.nan_to_num(companions['phase'],
                                                 nan=-99)
        bad = np.where((companions_phase_non_nan > 5) &
                       (companions_phase_non_nan < 101) &
                       (companions_phase_non_nan != -99))
        if self.verbose:
            print("Running the changing phase on companions")
            for ii in range(len(bad[0])):
                print('WARNING: changing phase {0} to 5'.format(companions_phase_non_nan[bad[0][ii]]))
            companions['phase'][bad] = 5
        # Get rid of the columns designation and and bad_system
        star_systemsPrime.remove_columns(['bad_system', 'designation', 'IFMR_it'])
        companions.remove_columns(['bad_system', 'IFMR_it'])
        # Following is commented out as I will need to rethink deletion.
        # if self.verbose:
        #.   print("{} non-single star systems".format(str(rejected_prims)) + 
        #.   " had to be deleted" +
        #.    " before IFMR application")
        #.    print("{} companions".format(str(rejected_companions)) +
        #.   " had to be deleted before" +
        #.   " IFMR was applied")
        # For testing purposes we can make rejected_system and rejected_comapnions
        # be returned too.
        return companions, star_systemsPrime

    def _remove_bad_systems(self, star_systems, companions=None):
        """
        Helper function to remove stars with masses outside the isochrone
        mass range from the cluster. These stars are identified by having 
        a Teff = 0, as set up by _make_star_systems_table_interp.
        If self.ifmr == None, then both high and low-mass bad systems are 
        removed. If self.ifmr != None, then we will save the high mass systems 
        since they will be plugged into an ifmr later.
        
        Unlike this function's equivalent in the ResolvedCluster, I directly use the minimum
        and maximum mass. This is since Teff NaN in the secondary stars may mean that
        the 
        """
        N_systems = len(star_systems)

        # Get rid of the bad ones
        # Convert nan_to_num to avoid errors on greater than, less than comparisons
        star_systems_mass_non_nan = np.nan_to_num(star_systems['mass'], nan=0)
        star_systems_phase_non_nan = np.nan_to_num(star_systems['phase'], nan=0)
        if self.ifmr == None:
            # Trim out bad systems; specifically, stars with masses outside those provided
            # by the model isochrone (except for compact objects, which is where I check
            # phases and allow stars with 101 + phase to not get filtered out.
            idx = np.where(((star_systems_mass_non_nan >= self.iso.min_mass) &
                           (star_systems_mass_non_nan <= self.iso.max_mass)) |
                           (star_systems_phase_non_nan >= 101))[0]
            # As in the old _remove_bad_systems we make sure that the for the
            # very massive stars, the log_g is applied
        else:
            # We want to trim out bad systems; specifically, stars with masses outside
            # those provided by the model isochrone (except for compact objects,
            # which is where I check phases and allow stars with 101 + phase to not get
            # filtered out.
            idx = np.where((star_systems_mass_non_nan >= self.iso.min_mass)|
                           (star_systems_phase_non_nan >= 101))[0]

        if len(idx) != N_systems and self.verbose:
            print( 'Found {0:d} stars out of mass range'.format(N_systems - len(idx)))        
        star_systems = star_systems[idx]
        # In the case IFMR!=None, we must designate stars with
        # mass above max_mass of isochrone to be inputted into the
        # IFMR.
        to_still_apply_IFMR = np.where(star_systems['mass'] > self.iso.max_mass)[0]
        if to_still_apply_IFMR and self.ifmr != None:
            star_systems['IFMR_it'][to_still_apply_IFMR] = True
        N_systems = len(star_systems)

        if companions!=None:
            # Clean up companion stuff (which we haven't handled yet)
            # Delete stars that are companions of primaries we deleted
            companions = companions[np.where(np.isin(star_systems['designation'], idx))[0]]
        
        return star_systems, companions

class ResolvedCluster(Cluster):
    """
    Cluster sub-class that produces a *resolved* stellar cluster.
    A table is output with the synthetic photometry and intrinsic 
    properties of the individual stars (or stellar systems, if 
    mutliplicity is used in the IMF object).
    If multiplicity is used, than a second table is produced that 
    contains the properties of the companion stars independent of their
    primary stars.
    Parameters
    -----------
    iso: isochrone object
        SPISEA isochrone object
    
    imf: imf object
        SPISEA IMF object
    cluster_mass: float
        Total initial mass of the cluster, in M_sun
    ifmr: ifmr object or None
        If ifmr object is defined, will create compact remnants
        produced by the cluster at the given isochrone age. Otherwise,
        no compact remnants are produced.
    seed: int
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output.
        Default None
    vebose: boolean
        True for verbose output.
    """
    def __init__(self, iso, imf, cluster_mass, ifmr=None, verbose=True,
                     seed=None):
        Cluster.__init__(self, iso, imf, cluster_mass, ifmr=ifmr, verbose=verbose,
                             seed=seed)

        # Provide a user warning is random seed is set
        if seed is not None:
            print('WARNING: random seed set to %i' % seed)

        t1 = time.time()
        ##### 
        # Sample the IMF to build up our cluster mass.
        #####
        mass, isMulti, compMass, sysMass = imf.generate_cluster(cluster_mass,
                                                                    seed=seed)

        # Figure out the filters we will make.
        self.filt_names = self.set_filter_names()
        self.cluster_mass = cluster_mass

        #####
        # Make isochrone interpolators
        #####
        interp_keys = ['Teff', 'L', 'logg', 'isWR', 'mass_current', 'phase'] + self.filt_names
        self.iso_interps = {}
        for ikey in interp_keys:
            self.iso_interps[ikey] = interpolate.interp1d(self.iso.points['mass'], self.iso.points[ikey],
                                                          kind='linear', bounds_error=False, fill_value=np.nan)
        
        ##### 
        # Make a table to contain all the information about each stellar system.
        #####
        star_systems = self._make_star_systems_table(mass, isMulti, sysMass)
        
        # Trim out bad systems; specifically, stars with masses outside those provided
        # by the model isochrone (except for compact objects).
        star_systems, compMass = self._remove_bad_systems(star_systems, compMass)

        ##### 
        # Make a table to contain all the information about companions.
        #####
        if self.imf.make_multiples:
            companions = self._make_companions_table(star_systems, compMass)
            
        #####
        # Save our arrays to the object
        #####
        self.star_systems = star_systems
        
        if self.imf.make_multiples:
            self.companions = companions

        return

    def set_filter_names(self):
        """
        Set filter column names
        """
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

        # Add columns for the Teff, L, logg, isWR, mass_current, phase for the primary stars.
        star_systems.add_column( Column(np.zeros(N_systems, dtype=float), name='Teff') )
        star_systems.add_column( Column(np.empty(N_systems, dtype=float), name='L') )
        star_systems.add_column( Column(np.empty(N_systems, dtype=float), name='logg') )
        star_systems.add_column( Column(np.empty(N_systems, dtype=float), name='isWR') )
        star_systems.add_column( Column(np.empty(N_systems, dtype=float), name='mass_current') )
        star_systems.add_column( Column(np.empty(N_systems, dtype=float), name='phase') )
        star_systems.add_column( Column(np.empty(N_systems, dtype=float), name='metallicity') )

        # Add the filter columns to the table. They are empty so far.
        # Keep track of the filter names in : filt_names
        for filt in self.filt_names:
            star_systems.add_column( Column(np.empty(N_systems, dtype=float), name=filt) )

        # Use our pre-built interpolators to fetch values from the isochrone for each star.
        star_systems['Teff'] = self.iso_interps['Teff'](star_systems['mass'])
        star_systems['L']    = self.iso_interps['L'](star_systems['mass'])
        star_systems['logg'] = self.iso_interps['logg'](star_systems['mass'])
        star_systems['isWR'] = np.round(self.iso_interps['isWR'](star_systems['mass']))
        star_systems['mass_current'] = self.iso_interps['mass_current'](star_systems['mass'])
        star_systems['phase'] = np.round(self.iso_interps['phase'](star_systems['mass']))
        star_systems['metallicity'] = np.ones(N_systems) * self.iso.metallicity

        # For a very small fraction of stars, the star phase falls on integers in-between
        # the ones we have definition for, as a result of the interpolation. For these
        # stars, round phase down to nearest defined phase (e.g., if phase is 71,
        # then round it down to 5, rather than up to 101).
        # Note: this only becomes relevant when the cluster is > 10**6 M-sun, this
        # effect is so small
        # Convert nan_to_num to avoid errors on greater than, less than comparisons
        star_systems_phase_non_nan = np.nan_to_num(star_systems['phase'], nan=-99)
        bad = np.where((star_systems_phase_non_nan > 5) & (star_systems_phase_non_nan < 101) &
                       (star_systems_phase_non_nan != 9) & (star_systems_phase_non_nan != -99))
        # Print warning, if desired
        verbose=False
        if verbose:
            for ii in range(len(bad[0])):
                print('WARNING: changing phase {0} to 5'.format(star_systems['phase'][bad[0][ii]]))
        star_systems['phase'][bad] = 5
        
        for filt in self.filt_names:
            star_systems[filt] = self.iso_interps[filt](star_systems['mass'])

        #####
        # Make Remnants
        #     Note: Some models already have WDs in them. If they do, then they shouldn't
        #     be handled by this code here (because their Teff > 0).
        # 
        # Remnants have flux = 0 in all bands if they are generated here.
        ##### 
        if self.ifmr != None:
            # Identify compact objects as those with Teff = 0 or with phase > 100.
            highest_mass_iso = self.iso.points['mass'].max()
            idx_rem = np.where((np.isnan(star_systems['Teff'])) & (star_systems['mass'] > highest_mass_iso))[0]
            
            # Calculate remnant mass and ID for compact objects; update remnant_id and
            # remnant_mass arrays accordingly
            if 'metallicity_array' in inspect.getfullargspec(self.ifmr.generate_death_mass).args:
                r_mass_tmp, r_id_tmp = self.ifmr.generate_death_mass(mass_array=star_systems['mass'][idx_rem],
                                                                     metallicity_array=star_systems['metallicity'][idx_rem])
            else:
                r_mass_tmp, r_id_tmp = self.ifmr.generate_death_mass(mass_array=star_systems['mass'][idx_rem])

            # Drop remnants where it is not relevant (e.g. not a compact object or
            # outside mass range IFMR is defined for)
            good = np.where(r_id_tmp > 0)
            idx_rem_good = idx_rem[good]

            star_systems['mass_current'][idx_rem_good] = r_mass_tmp[good]
            star_systems['phase'][idx_rem_good] = r_id_tmp[good]

            # Give remnants a magnitude of nan, so they can be filtered out later when calculating flux.
            for filt in self.filt_names:
                star_systems[filt][idx_rem_good] = np.full(len(idx_rem_good), np.nan)

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

        # Add columns for the Teff, L, logg, isWR mass_current, phase, and filters for the companion stars.
        companions.add_column( Column(np.zeros(N_comp_tot, dtype=float), name='mass') )
        companions.add_column( Column(np.zeros(N_comp_tot, dtype=float), name='Teff') )
        companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name='L') )
        companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name='logg') )
        companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name='isWR') )
        companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name='mass_current') )
        companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name='phase') )
        companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name='metallicity') )
        for filt in self.filt_names:
            companions.add_column( Column(np.empty(N_comp_tot, dtype=float), name=filt) )
            
        if isinstance(self.imf._multi_props, multiplicity.MultiplicityResolvedDK):
            companions.add_column( Column(np.zeros(N_comp_tot, dtype=float), name='log_a') )
            companions.add_column( Column(np.zeros(N_comp_tot, dtype=float), name='e') )
            companions.add_column( Column(np.zeros(N_comp_tot, dtype=float), name='i', description = 'degrees') )
            companions.add_column( Column(np.zeros(N_comp_tot, dtype=float), name='Omega') )
            companions.add_column( Column(np.zeros(N_comp_tot, dtype=float), name='omega') )
            
            for ii in range(len(companions)):
                companions['log_a'][ii] = self.imf._multi_props.log_semimajoraxis(star_systems['mass'][companions['system_idx'][ii]])
            companions['e'] = self.imf._multi_props.random_e(np.random.rand(N_comp_tot))
            companions['i'], companions['Omega'], companions['omega'] = \
            self.imf._multi_props.random_keplarian_parameters(np.random.rand(N_comp_tot),
                                                              np.random.rand(N_comp_tot),
                                                              np.random.rand(N_comp_tot))


        # Make an array that maps system index (ii), companion index (cc) to
        # the place in the 1D companions array.
        N_comp_max = N_companions.max()
        
        comp_index = np.zeros((N_systems, N_comp_max), dtype=int)
        kk = 0
        for ii in range(N_systems):
            for cc in range(N_companions[ii]):
                comp_index[ii][cc] = kk
                kk += 1

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
                companions['Teff'][cdx] = self.iso_interps['Teff'](comp_mass)
                companions['L'][cdx] = self.iso_interps['L'](comp_mass)
                companions['logg'][cdx] = self.iso_interps['logg'](comp_mass)
                companions['isWR'][cdx] = np.round(self.iso_interps['isWR'](comp_mass))
                companions['mass_current'] = self.iso_interps['mass_current'](companions['mass'])
                companions['phase'] = np.round(self.iso_interps['phase'](companions['mass']))
                companions['metallicity'] = np.ones(N_comp_tot)*self.iso.metallicity

                # For a very small fraction of stars, the star phase falls on integers in-between
                # the ones we have definition for, as a result of the interpolation. For these
                # stars, round phase down to nearest defined phase (e.g., if phase is 71,
                # then round it down to 5, rather than up to 101).
                # Convert nan_to_num to avoid errors on greater than, less than comparisons
                star_systems_phase_non_nan = np.nan_to_num(star_systems['phase'], nan=-99)
                bad = np.where((star_systems_phase_non_nan > 5) & (star_systems_phase_non_nan < 101) &
                               (star_systems_phase_non_nan != 9) & (star_systems_phase_non_nan != -99))
                # Print warning, if desired
                verbose=False
                if verbose:
                    for ii in range(len(bad[0])):
                        print('WARNING: changing phase {0} to 5'.format(companions['phase'][bad[0][ii]]))
                companions['phase'][bad] = 5

                for filt in self.filt_names:
                    # Magnitude of companion
                    companions[filt][cdx] = self.iso_interps[filt](comp_mass)

                    mag_s = star_systems[filt][idx]
                    mag_c = companions[filt][cdx]

                    # Add companion flux to system flux.
                    f1 = 10**(-mag_s / 2.5)
                    f2 = 10**(-mag_c / 2.5)

                    # For dark objects, turn the np.nan fluxes into zeros.
                    f1 = np.nan_to_num(f1)
                    f2 = np.nan_to_num(f2)

                    # If *both* objects are dark, then keep the magnitude
                    # as np.nan. Otherwise, add fluxes together
                    good = np.where( (f1 != 0) | (f2 != 0) )
                    bad = np.where( (f1 == 0) & (f2 == 0) )
                    
                    star_systems[filt][idx[good]] = -2.5 * np.log10(f1[good] + f2[good])
                    star_systems[filt][idx[bad]] = np.nan

        #####
        # Make Remnants with flux = 0 in all bands.
        ##### 
        if self.ifmr != None:
            # Identify compact objects as those with Teff = 0 or with masses above the max iso mass
            highest_mass_iso = self.iso.points['mass'].max()
            cdx_rem = np.where((companions['Teff'] == 0) &
                                (companions['mass'] > highest_mass_iso))[0]
            
            # Calculate remnant mass and ID for compact objects; update remnant_id and
            # remnant_mass arrays accordingly
            if 'metallicity_array' in inspect.getfullargspec(self.ifmr.generate_death_mass).args:
                r_mass_tmp, r_id_tmp = self.ifmr.generate_death_mass(mass_array=star_systems['mass'][cdx_rem],
                                                                     metallicity_array=star_systems['metallicity'][cdx_rem])
            else:
                r_mass_tmp, r_id_tmp = self.ifmr.generate_death_mass(mass_array=star_systems['mass'][cdx_rem])

            # Drop remnants where it is not relevant (e.g. not a compact object or
            # outside mass range IFMR is defined for)
            good = np.where(r_id_tmp > 0)
            cdx_rem_good = cdx_rem[good]

            companions['mass_current'][cdx_rem_good] = r_mass_tmp[good]
            companions['phase'][cdx_rem_good] = r_id_tmp[good]
                
            # Give remnants a magnitude of nan, so they can be filtered out later when calculating flux.
            for filt in self.filt_names:
                companions[filt][cdx_rem_good] = np.full(len(cdx_rem_good), np.nan)


        # Notify if we have a lot of bad ones.
        # Convert nan_to_num to avoid errors on greater than, less than comparisons
        companions_teff_non_nan = np.nan_to_num(companions['Teff'], nan=-99)
        idx = np.where(companions_teff_non_nan > 0)[0]
        if len(idx) != N_comp_tot and self.verbose:
            print( 'Found {0:d} companions out of stellar mass range'.format(N_comp_tot - len(idx)))

        # Double check that everything behaved properly.
        assert companions['mass'][idx].min() > 0

        return companions

    
    def _remove_bad_systems(self, star_systems, compMass):
        """
        Helper function to remove stars with masses outside the isochrone
        mass range from the cluster. These stars are identified by having 
        a Teff = 0, as set up by _make_star_systems_table_interp.
        If self.ifmr == None, then both high and low-mass bad systems are 
        removed. If self.ifmr != None, then we will save the high mass systems 
        since they will be plugged into an ifmr later.
        """
        N_systems = len(star_systems)

        # Get rid of the bad ones
        # Convert nan_to_num to avoid errors on greater than, less than comparisons
        star_systems_teff_non_nan = np.nan_to_num(star_systems['Teff'], nan=-99)
        star_systems_phase_non_nan = np.nan_to_num(star_systems['phase'], nan=-99)
        if self.ifmr == None:
            # Keep only those stars with Teff assigned.
            idx = np.where(star_systems_teff_non_nan > 0)[0]
        else:
            # Keep stars (with Teff) and any other compact objects (with phase info). 
            idx = np.where( (star_systems_teff_non_nan > 0) | (star_systems_phase_non_nan >= 0) )[0]

        if len(idx) != N_systems and self.verbose:
            print( 'Found {0:d} stars out of mass range'.format(N_systems - len(idx)))

        star_systems = star_systems[idx]
        N_systems = len(star_systems)

        if self.imf.make_multiples:
            # Clean up companion stuff (which we haven't handled yet)
            compMass = [compMass[ii] for ii in idx]
        
        return star_systems, compMass
    
class ResolvedClusterDiffRedden(ResolvedCluster):
    
    
    """
    Sub-class of ResolvedCluster that applies differential
    extinction to the synthetic photometry.

    Parameters
    -----------
    iso: isochrone object
        SPISEA isochrone object
    
    imf: imf object
        SPISEA IMF object

    cluster_mass: float
        Total initial mass of the cluster, in M_sun

    delta_AKs: float
        Amount of differential extinction to apply to synthetic photometry,
        in terms of magnitudes of extinction in the Ks filter. Specifically,
        delta_AKs defines the standard deviation of a Gaussian distribution 
        from which the delta_AKs values will be drawn from for each individual
        system.

    ifmr: ifmr object or None
        If ifmr object is defined, will create compact remnants
        produced by the cluster at the given isochrone age. Otherwise,
        no compact remnants are produced.

    seed: int
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output.
        Default None

    vebose: boolean
        True for verbose output.
    """
    def __init__(self, iso, imf, cluster_mass, deltaAKs,
                 ifmr=None, verbose=False, seed=None):

        
        ResolvedCluster.__init__(self, iso, imf, cluster_mass, ifmr=ifmr, verbose=verbose,
                                     seed=seed)

        # Set random seed, if desired
        if seed is not None:
            np.random.seed(seed=seed)

        # Extract the extinction law from the isochrone object
        redlaw_str = iso.points.meta['REDLAW']
        red_law = reddening.get_red_law(redlaw_str)
            
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
        # Aks +/- deltaAKs in each filter
        rand_red = np.random.randn(len(self.star_systems))

        for filt in self.filt_names:
            self.star_systems[filt] += rand_red * delta_red_filt[filt]

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
    
    
class UnresolvedCluster(Cluster):
    
    
    """
    Cluster sub-class that produces an *unresolved* stellar cluster.
    Output is a combined spectrum that is the sum of the individual 
    spectra of the cluster stars.

    Parameters
    -----------
    iso: isochrone object
        SPISEA isochrone object
    
    imf: imf object
        SPISEA IMF object

    cluster_mass: float
        Total initial mass of the cluster, in M_sun

    wave_range: 2-element array
        Define the minumum and maximum wavelengths of the final
        output spectrum, in Angstroms. Array should be [min_wave, max_wave]

    vebose: boolean
        True for verbose output.
    """
    def __init__(self, iso, imf, cluster_mass,
                 wave_range=[3000, 52000], verbose=False):
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
        self.wave_trim = self.spec_list_trim[0].wave
        
        t4 = time.time()
        print( 'Spec trimming took {0:f}s'.format(t4-t3))

        self.mass_tot = np.sum(sysMass[idx])
        print( 'Total cluster mass is {0:f} M_sun'.format(self.mass_tot))

        return
        
        
class Isochrone(object):
    
    
    """
    Base Isochrone class. 

    Parameters
    ----------
    logAge : float
        The age of the isochrone, in log(years)

    AKs : float
        The total extinction in Ks filter, in magnitudes

    distance : float
        The distance of the isochrone, in pc

    metallicity : float, optional
        The metallicity of the isochrone, in [M/H].
        Default is 0.

    evo_model: model evolution class, optional
        Set the stellar evolution model class. 
        Default is evolution.MISTv1().

    atm_func: model atmosphere function, optional
        Set the stellar atmosphere models for the stars. 
        Default is get_merged_atmosphere.

    wd_atm_func: white dwarf model atmosphere function, optional
        Set the stellar atmosphere models for the white dwafs. 
        Default is get_wd_atmosphere   

    mass_sampling : int, optional
        Sample the raw isochrone every `mass_sampling` steps. The default
        is mass_sampling = 0, which is the native isochrone mass sampling 
        of the evolution model.

    wave_range : list, optional
        length=2 list with the wavelength min/max of the final spectra.
        Units are Angstroms. Default is [3000, 52000].

    min_mass : float or None, optional
        If float, defines the minimum mass in the isochrone.
        Unit is solar masses. Default is None

    max_mass : float or None, optional
        If float, defines the maxmimum mass in the isochrone.
        Units is solar masses. Default is None.

    rebin : boolean, optional
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is False,
        which is often sufficient synthetic photometry in most cases.
    """
    
    
    def __init__(self, logAge, AKs, distance, metallicity=0.0,
                 evo_model=default_evo_model, atm_func=default_atm_func,
                 wd_atm_func = default_wd_atm_func,
                 red_law=default_red_law, mass_sampling=1,
                 wave_range=[3000, 52000], min_mass=None, max_mass=None,
                 rebin=True):


        t1 = time.time()
        
        c = constants

        # Assert that the wavelength ranges are within the limits of the
        # VEGA model (0.1 - 10 microns)
        try:
            assert wave_range[0] > 1000
            assert wave_range[1] < 100000
        except:
            print('Desired wavelength range invalid. Limit to 1000 - 10000 A')
            return
        
        # Get solar metallicity models for a population at a specific age.
        # Takes about 0.1 seconds.
        evol = evo_model.isochrone(age=10**logAge,
                                   metallicity=metallicity)

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

        # Give luminosity, temperature, mass, radius units (astropy units).
        L_all = 10**evol['logL'] * constants.L_sun # luminsoity in W
        T_all = 10**evol['logT'] * units.K
        R_all = np.sqrt(L_all / (4.0 * math.pi * c.sigma_sb * T_all**4))
        mass_all = evol['mass'] * units.Msun # masses in solar masses
        logg_all = evol['logg'] # in cgs
        mass_curr_all = evol['mass_current'] * units.Msun
        phase_all = evol['phase']
        isWR_all = evol['isWR']

        # Define the table that contains the "average" properties for each star.
        tab = Table([L_all, T_all, R_all, mass_all, logg_all, isWR_all, mass_curr_all, phase_all],
                    names=['L', 'Teff', 'R', 'mass', 'logg', 'isWR', 'mass_current', 'phase'])

        # Initialize output for stellar spectra
        self.spec_list = []

        # For each temperature extract the synthetic photometry.
        for ii in range(len(tab['Teff'])):
            # Loop is currently taking about 0.11 s per iteration
            gravity = float( logg_all[ii] )
            L = float( L_all[ii].cgs / (units.erg / units.s)) # in erg/s
            T = float( T_all[ii] / units.K)               # in Kelvin
            R = float( R_all[ii].to('pc') / units.pc)              # in pc
            phase = phase_all[ii]

            # Get the atmosphere model now. Wavelength is in Angstroms
            # This is the time-intensive call... everything else is negligable.
            # If source is a star, pull from star atmospheres. If it is a WD,
            # pull from WD atmospheres
            if phase == 101:
                star = wd_atm_func(temperature=T, gravity=gravity, metallicity=metallicity,
                                       verbose=False)
            else:
                star = atm_func(temperature=T, gravity=gravity, metallicity=metallicity,
                                    rebin=rebin)

            # Trim wavelength range down to JHKL range (0.5 - 5.2 microns)
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
        tab.meta['METAL_IN'] = evol.meta['metallicity_in']
        tab.meta['METAL_ACT'] = evol.meta['metallicity_act']
        tab.meta['WAVEMIN'] = wave_range[0]
        tab.meta['WAVEMAX'] = wave_range[1]

        self.points = tab

        t2 = time.time()
        print( 'Isochrone generation took {0:f} s.'.format(t2-t1))
        return

    def plot_HR_diagram(self, savefile=None):
        
        
        """
        Make a standard HR diagram for this isochrone.

        Parameters
        -----------
        savefile: path or None, optional
             Path to file plot too, if desired. 
             Default is None
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

        Parameters
        -----------
        savefile: path or None, optional
             Path to file plot too, if desired. 
             Default is None
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


class Isochrone_Binary(Isochrone):
    """
    Base Isochrone class.

    Parameters
    ----------
    logAge : float
        The age of the isochrone, in log(years)
    AKs : float
        The total extinction in Ks filter, in magnitudes
    distance : float
        The distance of the isochrone, in pc
    metallicity : float, optional
        The metallicity of the isochrone, in [M/H].
        Default is 0.
    evo_model: model evolution class, optional
        Set the stellar evolution model class.
        Default is evolution.MISTv1().
    atm_func: model atmosphere function, optional
        Set the stellar atmosphere models for the stars.
        Default is get_merged_atmosphere.
    wd_atm_func: white dwarf model atmosphere function, optional
        Set the stellar atmosphere models for the white dwafs.
        Default is get_wd_atmosphere
    mass_sampling : int, optional
        Sample the raw isochrone every `mass_sampling` steps. The default
        is mass_sampling = 0, which is the native isochrone mass sampling
        of the evolution model.
    wave_range : list, optional
        length=2 list with the wavelength min/max of the final spectra.
        Units are Angstroms. Default is [3000, 52000].
    min_mass : float or None, optional
        If float, defines the minimum mass in the isochrone.
        Unit is solar masses. Default is None
    max_mass : float or None, optional
        If float, defines the maxmimum mass in the isochrone.
        Units is solar masses. Default is None.
    rebin : boolean, optional
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is False,
        which is often sufficient synthetic photometry in most cases.
    Important Attributes:
    Primaries-- AstroPy Table containing the following columns:
    -----------------------------------------------------------
    mass = in units of Solar Masses. Contains the initial mass of the primary
    star
    log_a = the log10(separation of the star with the secondary in AU)
    mass_current = in units of solar masses. Contains the initial mass of
    the primary star
    L = Luminosity of the primary star in Watts
    Teff -- in Kelvin. Effective temperature of the primary star
    R -- in units of solar radii. Radius of the primary star.
    phase -- integer indicating whether a secondary star is a
    white dwarf (101), neutron star (102), black hole, already merged star,
    or is none of the previous (5)
    gravity -- log10(acceleration of gravity at the primary star's surface
    in m/s^2)
    isWR -- boolean indicating whether a primary star is a WR Star
    ------------------------------------------------------------
    singles-- AstroPy Table containing the following columns:
    -----------------------------------------------------------
    mass = in units of Solar Masses. Contains the initial mass
    of the single star
    mass_current = in units of solar masses. Contains the
    initial mass of the single star
    L = Luminosity of the single star in Watts
    Teff -- in Kelvin. Effective temperature of the single star
    R -- in units of solar radii. Radius of the single star.
    phase -- integer indicating whether a secondary star is a
    white dwarf (101), neutron star (102), black hole, already merged star,
    or is none of the previous (5)
    gravity -- log10(acceleration of gravity at the single
    star's surface in m/s^2)
    isWR -- boolean indicating whether a single star is a WR Star
    ------------------------------------------------------------
    secondaries -- AstroPy Table containing the following columns:
    -----------------------------------------------------------
    mass = in units of Solar Masses. Contains the initial mass
    of the secondary star
    mass_current = in units of solar masses. Contains the initial
    mass of the secondary star
    L = Luminosity of the secondary star in Watts
    Teff -- in Kelvin. Effective temperature of the single star
    R -- in units of solar radii. Radius of the secondary star.
    phase -- integer indicating whether a secondary star is a
    white dwarf (101), neutron star (102), black hole, already merged star,
    or is none of the previous (5)
    logg -- log10(acceleration of gravity at the
    secondary star's surface in m/s^2)
    isWR -- boolean indicating whether a secondary star is a WR Star
    merged -- whether the corresponding primary star (WHICH IS AT THE
    SAME INDEX in the primaries table as is the secondary star),
    has actually incorporated the secondary.
    If merged is true, the L, T_eff, gravity, mass_current will be set
    to np.nan along with any magnitudes.
    ------------------------------------------------------------
    """

    def __init__(self, logAge, AKs, distance, metallicity,
                 evo_model=evolution.BPASS(), atm_func=default_atm_func,
                 wd_atm_func=default_wd_atm_func, mass_sampling=1,
                 red_law=default_red_law,
                 wave_range=[3000, 52000], min_mass=None, max_mass=None,
                 filters=['ubv,U', 'ubv,V', 'ubv,B', 'ubv,R', 'ubv,I'],
                 rebin=True):
        t1=time.time()
        self.metallicity = metallicity
        self.logage = logAge
        # Accounting for the definition of metallicity of the
        # evolution object's
        # Isochrone function.
        # Changes by Ryota: make atm_func and wd_atm_func instance vars
        self.atm_func = atm_func
        self.wd_atm_func = wd_atm_func
        self.distance = distance
        self.wave_range = wave_range
        self.AKs = AKs
        self.red_law = red_law
        self.filters = filters
        self.filt_names_table = []
        self.verbose = True

        # Assert that the wavelength ranges are within the limits of the
        # VEGA model (0.1 - 10 microns)
        try:
            assert wave_range[0] > 1000
            assert wave_range[1] < 100000
        except AssertionError:
            print('Desired wavelength range invalid. Limit to 1000 - 10000 A')
            return
        # Get solar metallicity models for a population at a specific age.
        # Takes about 0.1 seconds.
        evol = evo_model.isochrone(age=10 ** logAge,
                                   metallicity=metallicity)

        # Eliminate cases where log g is less than 0
        idx = np.where(evol['logg'] > 0)
        evol = evol[idx]

        # Trim to desired mass range
        if min_mass:
            idx = np.where((evol['mass'] >= min_mass) & (evol['mass2'] >= min_mass))
            evol = evol[idx]
        if max_mass:
            idx = np.where((evol['mass'] <= max_mass) & (evol['mass2'] <= min_mass))
            evol = evol[idx]
        self.max_mass = np.max([np.max(evol['mass']), np.max(evol['mass2'])])
        self.min_mass = np.min([np.min(evol['mass']), np.min(evol['mass2'])])
        # Give luminosity, temperature, mass, radius units (astropy units).
        L_all = 10 ** evol['logL'] * constants.L_sun  # luminsoity in W
        T_all = 10 ** evol['logT'] * units.K
        # TO DO: Conditionals to make sure I am using valid values

        R_all = np.sqrt(L_all / (4.0 * math.pi * c.sigma_sb * T_all ** 4))
        # masses in solar masses
        mass_all = evol['mass'] * units.Msun
        logg_all = evol['logg']
        mass_curr_all = evol['mass_current'] * units.Msun
        # We will keep track of the phase of the BPASS primary/
        # single star system
        phase_all = evol['phase']
        # We will keep track of the phase of the BPASS secondary
        phase_all2 = evol['phase2']
        # We will keep track of whether the star is WR or not
        isWR_all = evol['isWR']
        isWR2 = evol['isWR2']
        # Actual temperature of secondary in Kelvins
        Tef2 = (10 ** evol['log(T2)']) * units.K
        R2 = (10 ** evol['log(R2)']) * constants.R_sun
        L2 = (10 ** evol['log(L2)']) * constants.L_sun
        singles = Table([mass_all, L_all, T_all, R_all, logg_all, isWR_all,
                        mass_curr_all, phase_all, evol['single'], evol['source']],
                        names=['mass', 'L', 'Teff', 'R', 'logg',
                               'isWR', 'mass_current', 'phase', 'single', 'source'])
        # Also have inserted conversion factor to deal with the units of
        # log(a)
        primaries = Table([mass_all,
                           L_all, T_all, R_all, logg_all, isWR_all,
                           mass_curr_all, phase_all, evol['single'], evol['source']],
                          names=['mass', 'L', 'Teff', 'R', 'logg',
                                 'isWR', 'mass_current', 'phase',  'single', 'source'])
        # Note that we keep information regarding whether
        # a star corresponds to a merger.
        # This will help us decide whether we should not account
        # for the L2, Tef_2 during atmosphere creation
        # Yes I do unit conversions for log(a) column
        # as that is in log(a/R_sun)
        secondaries = Table([evol['mass2'] * units.Msun,
                            evol['log(a)'] + np.log10(constants.R_sun/constants.au),
                            L2, Tef2, R2,
                            evol['logg2'], isWR2,
                            evol['mass_current2'] * units.Msun,
                            phase_all2,
                            evol['single'],
                            evol['mergered?'], evol['source']],
                            names=['mass', 'log_a', 'L', 'Teff',
                                   'R', 'logg', 'isWR', 'mass_current',
                                   'phase', 'single', 'merged', 'source'])
        # Make sure that we are only looking at stars with companions when
        # examining the secondaries and primaries.
        secondaries = secondaries[np.where(~secondaries['single'])[0]]
        secondaries.remove_column('single')
        singles = singles[np.where(singles['single'])[0]]
        singles.remove_column('single')
        primaries = primaries[np.where(~primaries['single'])[0]]
        primaries.remove_column('single')
        # Trim down the table by selecting every Nth point where
        # N = mass sampling factor.
        singles = singles[::mass_sampling]
        primaries = primaries[::mass_sampling]
        secondaries = secondaries[::mass_sampling]
        # Inserting before I forget
        # I try to make sure that we have a queue
        # of atm_function results to process
        # If we have null values
        self.spec_list_si = [None for x in range(len(singles))]  # For single Stars
        self.spec_list2_pri = [None for x in range(len(primaries))]  # For primary stars
        self.spec_list3_sec = [None for x in range(len(secondaries))]  # For secondary stars
        # Turns into an attribute since we will access this in another function
        self.pairings2 = {"Singles": self.spec_list_si,
                          "Primaries": self.spec_list2_pri,
                          "Secondaries": self.spec_list3_sec}
        pairings = {"Singles": singles, "Primaries": primaries,
                    "Secondaries": secondaries}
        codes = {"Singles": 0, "Primaries": 1,
                    "Secondaries": 2}
        self.codes2 = {0: singles, 1: primaries,
                       2: secondaries}
        self.pairings = pairings
        self.singles = singles
        self.primaries = primaries
        self.secondaries = secondaries
        # For each temperature extract the synthetic photometry.
        for x in pairings:
            tab = pairings[x]
            atm_list = self.pairings2[x]
            # Workaround for a glitch encountered with units not showing up.
            # may need to come back and get rid of it since it looks silly.
            R_all = tab['R'] * units.m/units.m
            gravity_table = tab['logg']
            # a little issue with the compact remnant primaries from
            # the secondary star
            if x == 'Secondaries':
                merged = np.where(tab['merged'])
                tab['mass_current'][merged] = np.nan
                tab['Teff'][merged] = np.nan
                tab['R'][merged] = np.nan
                tab['logg'][merged] = np.nan
                tab['isWR'][merged] = False
                tab['phase'][merged] = -99
                                 
            cond = np.where(np.isfinite(tab['logg']) & (tab['logg'] != 0.0) &
                            np.isfinite(tab['L']) & (tab['L'] > 0.0) &
                            np.isfinite(tab['Teff']) & (tab['Teff'] > 0.0) &
                            np.isfinite(tab['R']) & (tab['R']> 0.0) &
                            (tab['phase'] <= 101) & (tab['phase'] != -99))[0]
            
            wrapper = BypasserContainer()
            wrapper.wave_range = wave_range
            wrapper.distance = distance
            wrapper.lis = atm_list
            vectorized_atm_maker = np.vectorize(self.atm_generator_to_vectorize)
            vectorized_atm_maker(cond, wrapper, wrapper, metallicity, atm_func,
                                 wd_atm_func, red_law, AKs, rebin, codes[x])
            bad_starcond = np.where(~ (np.isfinite(tab['logg']) &
                                       (tab['logg'] != 0.0) &
                                       np.isfinite(tab['L']) &
                                       (tab['L'] > 0.0) &
                                       np.isfinite(tab['Teff']) &
                                       (tab['Teff'] > 0.0) &
                                       np.isfinite(tab['R']) & (tab['R'] > 0.0) &
                                       (tab['phase'] <= 101) &
                                       (tab['phase'] != -99)))[0]
            tab['Teff'][bad_starcond] = np.nan
            tab['L'][bad_starcond] = np.nan
            tab['R'][bad_starcond] = np.nan
            tab['logg'][bad_starcond] = np.nan
        # I hope to change the next few lines to this as this will make
        # more thorough use of numpy and decrease the number of for-loop
        # iterations that are used.
                      
        self.singles = singles
        self.primaries = primaries
        self.secondaries = secondaries
        # Append all the meta data to the summary table.
        for tab in (singles, primaries, secondaries):
            tab.meta['REDLAW'] = red_law.name
            tab.meta['ATMFUNC'] = atm_func.__name__
            tab.meta['EVOMODEL'] = 'BPASS v2.2'
            tab.meta['LOGAGE'] = logAge
            tab.meta['AKS'] = AKs
            tab.meta['DISTANCE'] = distance
            tab.meta['METAL_IN'] = evol.meta['metallicity_in']
            tab.meta['METAL_ACT'] = evol.meta['metallicity_act']
            tab.meta['WAVEMIN'] = wave_range[0]
            tab.meta['WAVEMAX'] = wave_range[1]
        self.make_photometry()
        t2 = time.time()
        print('Isochrone generation took {0:f} s.'.format(t2-t1))
        return

    def atm_generator_to_vectorize(self, c_ind, w_wrapper, list_wrapper, met,
                                   atm_func, wd_atm_func, red_law, AKs, rebin,
                                   code):
        """
        - Input -
        c_ind: 
        w_wrapper:
        list_wrapper:
        met: 
        rebin: 
        atm_func: atmosphere function inputted into the constructor of
        Isochrone_Binary
        red_law: Extinction law object inputted into constructor of isochrone
        binary
        AKs: extinction at the Ks band inputted 
        """
        tab = self.codes2[code]
        wave_range = w_wrapper.wave_range
        distance = w_wrapper.distance
        atm_list = list_wrapper.lis
        R_all = tab['R'] * units.m / units.m
        R = float( R_all[c_ind].to('pc') / units.pc) 

        if (tab[c_ind]['phase'] < 101):
            star =  atm_func(temperature=tab['Teff'][c_ind],
                             gravity=tab['logg'][c_ind],
                             metallicity=met,
                             rebin=rebin)
        else:
            star =  wd_atm_func(temperature=tab['Teff'][c_ind],
                                gravity=tab['logg'][c_ind],
                                metallicity=met,
                                verbose=False)
        # Trim wavelength range down to
        # JHKL range (0.5 - 5.2 microns)
        star = spectrum.trimSpectrum(star, wave_range[0],
                                             wave_range[1])
        # Convert into flux observed at Earth (unreddened)
        R = float(R_all[c_ind].to("pc") / units.pc)
        star *= (R / distance) ** 2 # in erg s^-1 cm^-2 A^-1
        # Redden the spectrum. This doesn't take much time at all.
        red = red_law.reddening(AKs).resample(star.wave)
        star *= red
        # Save the final spectrum to our spec_list for later use.
        atm_list[c_ind] = star
        return

    def plot_HR_diagram(self, savefile=None):
        
        
        """
        Make a standard HR diagram for this isochrone.
        Took the code from the Isochrone base class and
        adapted it to work on the isochrone binary object
        Parameters
        -----------
        savefile: path or None, optional
             Path to file plot too, if desired. 
             Default is None
        """
        plt.clf()
        plt.loglog(self.singles['Teff'], self.singles['L'],
                   color='black', linestyle='solid', marker='+')
        plt.loglog(self.primaries['Teff'], self.primaries['L'],
                   color='black', linestyle='solid', marker='+')
        plt.loglog(self.secondaries['Teff'], self.secondaries['L'],
                   color='black', linestyle='solid', marker='+')
        plt.gca().invert_xaxis()
        plt.xlabel(r'T$_{\mathrm{eff}}$ (K)')
        plt.ylabel('Luminosity (erg / s)')
        
        fmt_title = 'logAge={0:.2f}, d={1:.2f} kpc, AKs={2:.2f}'
        plt.title(fmt_title.format(self.singles.meta['LOGAGE'],
                                  self.singles.meta['DISTANCE']/1e3,
                                  self.singles.meta['AKS']))

        if savefile != None:
            plt.savefig(savefile)
        
        return

    def plot_mass_luminosity(self, savefile=None):
        
        
        """
        Make a standard mass-luminosity relation plot for this isochrone.
        Took the code from the Isochrone base class and
        adapted it to work on the isochrone binary object
        Parameters
        -----------
        savefile: path or None, optional
             Path to file plot too, if desired. 
             Default is None
        """
        plt.clf()
        plt.loglog(self.singles['mass'], self.singles['L'], 'k.')
        plt.loglog(self.primaries['mass'], self.primaries['L'], 'k.')
        plt.loglog(self.secondaries['mass'], self.secondaries['L'], 'k.')
        plt.xlabel(r'Mass (M$_\odot$)')
        plt.ylabel('Luminosity (erg / s)')
        
        fmt_title = 'logAge={0:.2f}, d={1:.2f} kpc, AKs={2:.2f}'
        plt.title(fmt_title.format(self.singles.meta['LOGAGE'],
                                  self.singles.meta['DISTANCE']/1e3,
                                  self.singles.meta['AKS']))
        
        if savefile != None:
            plt.savefig(savefile)
        
        return

    def make_photometry(self, rebin=True, vega=vega):
        """
        Make synthetic photometry for the specified filters. This function
        udpates the self.points table to include new columns with the
        photometry.
        """
        startTime = time.time()

        meta = self.singles.meta

        print('Making photometry for isochrone: log(t) = %.2f  AKs = %.2f  dist = %d' %(meta['LOGAGE'], meta['AKS'], meta['DISTANCE']))
        print('     Starting at: ', datetime.datetime.now(),
              '  Usually takes ~5 minutes')
        # npoints = len(self.points)
        verbose_fmt = 'M = {0:7.3f} Msun  T = {1:5.0f} K  m_{2:s} = {3:4.2f}'

        # Loop through the filters, get filter info, make photometry for
        # all stars in this filter.
        for ii in self.filters:
            prt_fmt = 'Starting filter: {0:s}   Elapsed time: {1:.2f} seconds'
            print(prt_fmt.format(ii, time.time() - startTime))
            filt = get_filter_info(ii, rebin=rebin, vega=vega)
            filt_name = get_filter_col_name(ii)
            # Make the column to hold magnitudes
            # in this filter. Add to points table.
            col_name = 'm_' + filt_name
            self.filt_names_table.append(col_name)
            mag_col = Column(np.zeros(len(self.singles), dtype=float),
                             name=col_name)
            self.singles.add_column(mag_col)
            mag_col = Column(np.zeros(len(self.secondaries), dtype=float),
                             name=col_name)
            self.secondaries.add_column(mag_col)
            mag_col = Column(np.zeros(len(self.primaries), dtype=float),
                             name=col_name)
            self.primaries.add_column(mag_col)
            # Loop through each star in the isochrone and
            # do the filter integration
            print('Starting synthetic photometry')
            for x in self.pairings:
                listofStars = self.pairings2[x]
                table = self.pairings[x]
                length_of_list = len(listofStars)
                for ss in range(length_of_list):
                    star = listofStars[ss]
                    # Nada used to be a fill in value.
                    # I thought that None would be more intuitive
                    # So we will be using None
                    if star:
                        # These are already extincted, observed spectra.
                        star_mag = mag_in_filter(star, filt)
                    else:
                        # We want nan to stand in place
                        # for stars that do not have valid T
                        # or radius or L.
                        star_mag = np.nan
                    table[col_name][ss] = star_mag
                    if (self.verbose and (ss % 100) == 0):
                        print(verbose_fmt.format(table['mass'][ss],
                                                 table['Teff'][ss],
                                                 filt_name, star_mag))
        endTime = time.time()
        print('      Time taken: {0:.2f} seconds'.format(endTime - startTime))
        return


class IsochronePhot(Isochrone):
    """
    Make an isochrone with synthetic photometry in various filters. 
    Load from file if possible. 

    Parameters
    ----------
    logAge : float
        The age of the isochrone, in log(years)

    AKs : float
        The total extinction in Ks filter, in magnitudes

    distance : float
        The distance of the isochrone, in pc

    metallicity : float, optional
        The metallicity of the isochrone, in [Z/Z_solar].
        Default is 0.

    evo_model: model evolution class, optional
        Set the stellar evolution model class. 
        Default is evolution.MISTv1().

    atm_func: model atmosphere function, optional
        Set the stellar atmosphere models for the stars. 
        Default is atmospheres.get_merged_atmosphere.

    wd_atm_func: white dwarf model atmosphere function, optional
        Set the stellar atmosphere models for the white dwafs. 
        Default is atmospheres.get_wd_atmosphere   

    red_law : reddening law object, optional
        Define the reddening law for the synthetic photometry.
        Default is reddening.RedLawNishiyama09().

    iso_dir : path, optional
         Path to isochrone directory. Code will check isochrone
         directory to see if isochrone file already exists; if it 
         does, it will just read the isochrone. If the isochrone 
         file doesn't exist, then save isochrone to the isochrone
         directory.

    mass_sampling : int, optional
        Sample the raw isochrone every `mass_sampling` steps. The default
        is mass_sampling = 0, which is the native isochrone mass sampling 
        of the evolution model.

    wave_range : list, optional
        length=2 list with the wavelength min/max of the final spectra.
        Units are Angstroms. Default is [3000, 52000].

    min_mass : float or None, optional
        If float, defines the minimum mass in the isochrone.
        Unit is solar masses. Default is None

    max_mass : float or None, optional
        If float, defines the maxmimum mass in the isochrone.
        Units is solar masses. Default is None.

    rebin : boolean, optional
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is True,
        which is often sufficient synthetic photometry in most cases.

    recomp : boolean, optional
        If true, recalculate the isochrone photometry even if 
        the savefile exists

    filters : array of strings, optional
        Define what filters the synthetic photometry
        will be calculated for, via the filter string 
        identifier. 
    """
    def __init__(self, logAge, AKs, distance,
                 metallicity=0.0,
                 evo_model=default_evo_model, atm_func=default_atm_func,
                 wd_atm_func = default_wd_atm_func,
                 wave_range=[3000, 52000],
                 red_law=default_red_law, mass_sampling=1, iso_dir='./',
                 min_mass=None, max_mass=None, rebin=True, recomp=False,
                 filters=['ubv,U', 'ubv,B', 'ubv,V',
                          'ubv,R', 'ubv,I']):

        self.metallicity = metallicity

        # Make the iso_dir, if it doesn't already exist
        if not os.path.exists(iso_dir):
            os.mkdir(iso_dir)

        # Make and input/output file name for the stored isochrone photometry.
        # For solar metallicity case, allow for legacy isochrones (which didn't have
        # metallicity tag since they were all solar metallicity) to be read
        # properly
        if metallicity == 0.0:
            save_file_fmt = '{0}/iso_{1:.2f}_{2:4.2f}_{3:4s}_p00.fits'
            self.save_file = save_file_fmt.format(iso_dir, logAge, AKs, str(distance).zfill(5))

            save_file_legacy = '{0}/iso_{1:.2f}_{2:4.2f}_{3:4s}.fits'
            self.save_file_legacy = save_file_legacy.format(iso_dir, logAge, AKs, str(distance).zfill(5))
        else:
            # Set metallicity flag
            if metallicity < 0:
                metal_pre = 'm'
            else:
                metal_pre = 'p'
            metal_flag = int(abs(metallicity)*10)
  
                                                                                     
            save_file_fmt = '{0}/iso_{1:.2f}_{2:4.2f}_{3:4s}_{4}{5:2s}.fits'
            self.save_file = save_file_fmt.format(iso_dir, logAge, AKs, str(distance).zfill(5), metal_pre, str(metal_flag).zfill(2))
            self.save_file_legacy = save_file_fmt.format(iso_dir, logAge, AKs, str(distance).zfill(5), metal_pre, str(metal_flag).zfill(2))
            
        # Expected filters
        self.filters = filters

        # Recalculate isochrone if save_file doesn't exist or recomp == True
        file_exists = self.check_save_file(evo_model, atm_func, red_law)

        if (not file_exists) | (recomp==True):
            self.recalc = True
            Isochrone.__init__(self, logAge, AKs, distance,
                               metallicity=metallicity,
                               evo_model=evo_model, atm_func=atm_func,
                               wd_atm_func=wd_atm_func,
                               wave_range=wave_range,
                               red_law=red_law, mass_sampling=mass_sampling,
                               min_mass=min_mass, max_mass=max_mass, rebin=rebin)
            self.verbose = True
            
            # Make photometry
            self.make_photometry(rebin=rebin, vega=vega)
        else:
            self.recalc = False
            try:
                self.points = Table.read(self.save_file)
            except:
                self.points = Table.read(self.save_file_legacy)
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
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                self.points.write(self.save_file, overwrite=True)

        return

    def check_save_file(self, evo_model, atm_func, red_law):
        """
        Check to see if save_file exists, as saved by the save_file 
        and save_file_legacy objects. If the filename exists, check the 
        meta-data as well.

        returns a boolean: True is file exists, false otherwise
        """
        out_bool = False
        
        if os.path.exists(self.save_file) | os.path.exists(self.save_file_legacy):
            try:
                tmp = Table.read(self.save_file)
            except:
                tmp = Table.read(self.save_file_legacy)
            
        
            # See if the meta-data matches: evo model, atm_func, redlaw
            if ( (tmp.meta['EVOMODEL'] == type(evo_model).__name__) &
                (tmp.meta['ATMFUNC'] == atm_func.__name__) &
                 (tmp.meta['REDLAW'] == red_law.name) ):
                out_bool = True
            
        return out_bool

    def plot_CMD(self, mag1, mag2, savefile=None):
        """
        Make a CMD with mag1 vs mag1 - mag2

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
# NOTE: THIS CLASS IS DEPRECATED, DO NOT USE!
#===================================================#
class iso_table(object):
    def __init__(self, logAge, distance, evo_model=default_evo_model,
                 atm_func=default_atm_func, mass_sampling=1,
                 min_mass=None, max_mass=None, wave_range=[5000, 52000],
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
        evo_model : SPISEA evolution object
            Stellar evolution models used
        atm_func: SPISEA atmosphere object
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

            

            # Get the atmosphere model now. Wavelength is in Angstroms
            # This is the time-intensive call... everything else is negligable.
            star = atm_func(temperature=T, gravity=gravity)
            
            # Trim wavelength range down to JHKL range (0.5 - 5.2 microns)
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
            
        extinction_law: SPISEA extinction object
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
    Vega spectrum. Input name is the SPISEA
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

    elif name.startswith('naco'):
        filt = filters.get_naco_filt(filterName)

    elif name.startswith('ubv'):
        filt = filters.get_ubv_filt(filterName)

    elif name.startswith('ukirt'):
        filt = filters.get_ukirt_filt(filterName)
        
    elif name.startswith('keck_osiris'):
        filt = filters.get_keck_osiris_filt(filterName)

    elif name.startswith('ztf'):
        filt = filters.get_ztf_filt(filterName)

    elif name.startswith('gaia'):
        version = tmp[1]
        filt = filters.get_gaia_filt(version, filterName)

    elif name.startswith('hawki'):
        filt = filters.get_hawki_filt(filterName)
        
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
    appropriate SPISEA obs_string
    """
    # How we deal with obs_string is slightly different depending
    # if it is an hst filter (and thus pysynphot syntax) or our
    # own defined filters
    tmp = obs_str.split(',')

    if len(tmp) == 3:
        # Catch Gaia filter cases. Otherwise, it is HST filter
        if 'dr2_rev' in tmp:
            filt_name = 'gaiaDR2_{0}'.format(tmp[-1])
        else:
            filt_name = 'hst_{0}'.format(tmp[-1])
    else:
        filt_name = '{0}_{1}'.format(tmp[0], tmp[1])
        
    return filt_name

def get_obs_str(col):
    """
    Helper function to get the associated SPISEA obs_str given
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
                 'jwst_F070W': 'jwst,F070W',
                 'jwst_F115W': 'jwst,F115W',
                 'jwst_F140M': 'jwst,F140M',
                 'jwst_F150W': 'jwst,F150W',
                 'jwst_F150W2': 'jwst,F150W2',
                 'jwst_F162M': 'jwst,F162M',
                 'jwst_F182M': 'jwst,F182M',
                 'jwst_F187N': 'jwst,F187N',
                 'jwst_F200W': 'jwst,F200W',
                 'jwst_F210M': 'jwst,F210M',
                 'jwst_F250M': 'jwst,F250M', 
                 'jwst_F277W': 'jwst,F277W',
                 'jwst_F300M': 'jwst,F300M',
                 'jwst_F322W2': 'jwst,F322W2',
                 'jwst_F335M': 'jwst,F335M',
                 'jwst_F356W': 'jwst,F356W',
                 'jwst_F360M': 'jwst,F360M',
                 'jwst_F405N': 'jwst,F405N',
                 'jwst_F410M': 'jwst,F410M',
                 'jwst_F430M': 'jwst,F430M',
                 'jwst_F444W': 'jwst,F444W',
                 'jwst_F440W': 'jwst,F440W',
                 'jwst_F460M': 'jwst,F460M',
                 'jwst_F470N': 'jwst,F470N',
                 'jwst_F480M': 'jwst,F480M',
                 'nirc2_J': 'nirc2,J', 'nirc2_H': 'nirc2,H', 'nirc2_Kp': 'nirc2,Kp', 'nirc2_K': 'nirc2,K',
                 'nirc2_Lp': 'nirc2,Lp', 'nirc2_Ms': 'nirc2,Ms', 'nirc2_Hcont': 'nirc2,Hcont',
                 'nirc2_FeII': 'nirc2,FeII', 'nirc2_Brgamma': 'nirc2,Brgamma',
                 '2mass_J': '2mass,J', '2mass_H': '2mass,H', '2mass_Ks': '2mass,Ks',
                 'ubv_U':'ubv,U', 'ubv_B':'ubv,B', 'ubv_V':'ubv,V', 'ubv_R':'ubv,R',
                 'ubv_I':'ubv,I', 
                 'jg_J': 'jg,J', 'jg_H': 'jg,H', 'jg_K': 'jg,K',
                 'nirc1_K':'nirc1,K', 'nirc1_H':'nirc1,H',
                 'naco_J':'naco,J', 'naco_H':'naco,H', 'naco_Ks':'naco,Ks',
                 'ukirt_J':'ukirt,J', 'ukirt_H':'ukirt,H', 'ukirt_K':'ukirt,K',
                 'ctio_osiris_H': 'ctio_osiris,H', 'ctio_osiris_K': 'ctio_osiris,K',
                 'ztf_g':'ztf,g', 'ztf_r':'ztf,r', 'ztf_i':'ztf,i',
                 'gaiaDR2_G': 'gaia,dr2_rev,G', 'gaiaDR2_Gbp':'gaia,dr2_rev,Gbp',
                 'gaiaDR2_Grp':'gaia,dr2_rev,Grp',
                 'hawki_J': 'hawki,J',
                 'hawki_H': 'hawki,H',
                 'hawki_Ks': 'hawki,Ks'}

    obs_str = filt_list[name]
        
    return obs_str

def rebin_spec(wave, specin, wavnew):
    """
    Helper function to rebin spectra, from Jessica Lu's post
    on Astrobetter:
    https://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/
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
 
    evo_models: SPISEA evolution object
        Which evolution models to use

    atm_models: SPISEA atmospheres object
        Which atmosphere models to use

    redlaw: SPISEA reddening object
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
    _out.write('SPISEA parameters used to generate isochrone grid:\n')
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

    # Model mass has to be within 10% of the desired mass
    if (dm[mdx] / theMass) > 0.1:
        return None
    else:
        return mdx

    
def get_evo_model_by_string(evo_model_string):
    return getattr(evolution, evo_model_string)


def calc_ab_vega_filter_conversion(filt_str):
    """
    Function to calculate the conversion between
    AB and Vega magnitudes for a given filter:
    m_AB - m_vega

    Note: this conversion is just the vega magnitude in 
    AB system

    Parameters:
    -----------
    filt_str: string
        Filter identification string
    """
    # Get filter info
    filt = get_filter_info(filt_str)

    # Let's convert everything into frequency space
    c = 2.997*10**18 # A / s
    vega_wave = vega.wave
    vega_mu = c / vega_wave
    vega_flux_mu = vega.flux * (vega_wave **2 / c)

    filt_wave = filt.wave
    filt_mu = c / filt_wave
    s_filt = filt.throughput
    
    # Interpolate the filter function, determine what the
    # filter function is at the exact sampling of the
    # vega spectrum (in freq space)
    filt_interp = scipy.interpolate.interp1d(filt_mu, s_filt, kind='linear', bounds_error=False,
                                                 fill_value=0)
    s_interp = filt_interp(vega_mu)
    
    # Now for the m_ab calculation
    mu_diff = np.diff(vega_mu)
    numerator = np.sum(vega_flux_mu[:-1] * s_interp[:-1] * mu_diff)
    denominator = np.sum(s_interp[:-1] * mu_diff)

    vega_mag_ab = -2.5 * np.log10(numerator / denominator) - 48.6

    print('For {0}, m_ab - m_vega = {1}'.format(filt_str, vega_mag_ab))

    #--Same calculation, in lambda space. Less accurate for some reason---#
    # Interpolate the filter function to be the exact same sampling as the
    # vega spectrum
    #c = 3*10**18 # A / s
    #filt_interp = scipy.interpolate.interp1d(filt.wave, filt.throughput, kind='linear', bounds_error=False,
    #                                             fill_value=0)
    #s_interp = filt_interp(vega.wave)

    # Calculate the numerator 
    #diff = np.diff(vega.wave)
    #numerator2 = np.sum((vega.wave[:-1]**2. / c) * vega.flux[:-1] * s_interp[:-1] * diff)
    
    # Now we need to intergrate the filter response for the denominator
    #denominator2 = np.sum(s_interp[:-1] * diff)

    # Calculate vega AB magnitude. This is the conversion
    #vega_mag_ab2 = -2.5 * np.log10(numerator2 / denominator2) - 48.6

    return vega_mag_ab

