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
from astropy import constants, units
from astropy.table import Table, Column, MaskedColumn
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
import inspect
import astropy.modeling

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

class Interpolator(object):
    def __init__(self, xp, yp):
        """Wrapper for np.interp to allow for pickling in multiprocessing.

        Parameters
        ----------
        xp: array-like
            x data points for interpolation
        yp: array-like
            y data points for interpolation
        """
        self.xp = xp
        self.yp = yp

    def __call__(self, x):
        return np.interp(x, self.xp, self.yp, left=np.nan, right=np.nan)

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
        self.rng = np.random.default_rng(self.seed)

        return

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
    def __init__(self, iso, imf, cluster_mass, ifmr=None, verbose=False,
                     seed=None):
        Cluster.__init__(self, iso, imf, cluster_mass, ifmr=ifmr, verbose=verbose,
                             seed=seed)
        # Provide a user warning is random seed is set
        if seed is not None and verbose:
            print('WARNING: random seed set to %i' % seed)

        #####
        # Sample the IMF to build up our cluster mass.
        #####
        # start0 = time.time()
        mass, isMulti, compMass, sysMass = imf.generate_cluster(cluster_mass)
        # end0 = time.time()
        # print('IMF sampling took {0:f} s.'.format(end0 - start0))

        # Figure out the filters we will make.
        self.filt_names = self.set_filter_names()
        self.cluster_mass = cluster_mass

        #####
        # Make isochrone interpolators
        #####
        interp_keys = ['Teff', 'L', 'logg', 'isWR', 'mass_current', 'phase'] + self.filt_names
        self.iso_interps = {}
        for ikey in interp_keys:
            # self.iso_interps[ikey] = interpolate.interp1d(self.iso.points['mass'], self.iso.points[ikey],
            #                                               kind='linear', bounds_error=False, fill_value=np.nan)
            self.iso_interps[ikey] = Interpolator(self.iso.points['mass'], self.iso.points[ikey])

        #####
        # Make a table to contain all the information about each stellar system.
        #####
        # start1 = time.time()
        star_systems = self._make_star_systems_table(mass, isMulti, sysMass)
        # end1 = time.time()
        # print('Star systems table took {0:f} s.'.format(end1 - start1))

        # Trim out bad systems; specifically, stars with masses outside those provided
        # by the model isochrone (except for compact objects).
        # start2 = time.time()
        star_systems, compMass = self._remove_bad_systems(star_systems, compMass)
        # end2 = time.time()
        # print('Bad system removal took {0:f} s.'.format(end2-start2))

        #####
        # Make a table to contain all the information about companions.
        #####
        if self.imf.make_multiples:
            # start3 = time.time()
            star_systems, companions = self._make_companions_table_new(star_systems, compMass)
            # end3 = time.time()
            # print('Companion table new took {0:f} s.'.format(end3 - start3))
            self.companions = companions

            # compMass = [
            #     [value for value, mask in zip(row, row_mask) if not mask]
            #     for row, row_mask in zip(compMass.data, compMass.mask)
            # ]
            # start3 = time.time()
            # star_systems, companions = self._make_companions_table(star_systems, compMass)
            # end3 = time.time()
            # print('Companion table took {0:f} s.'.format(end3-start3))
            # self.companions = companions

        #####
        # Save our arrays to the object
        #####
        self.star_systems = star_systems

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

        # Use our pre-built interpolators to fetch values from the isochrone for each star.
        for key in ['Teff', 'L', 'logg', 'mass_current']:
            star_systems.add_column(Column(self.iso_interps[key](star_systems['mass']), name=key))

        # Treat out-of-range mass as isWR=True
        star_systems.add_column(Column(~(self.iso_interps['isWR'](star_systems['mass']) < 0.5), name='isWR'))
        star_systems.add_column(Column(np.round(self.iso_interps['phase'](star_systems['mass'])), name='phase'))

        star_systems['metallicity'] = np.ones(N_systems) * self.iso.metallicity

        # Add the filter columns to the table. They are empty so far.
        # Keep track of the filter names in : filt_names
        for filt in self.filt_names:
            star_systems.add_column(Column(self.iso_interps[filt](star_systems['mass']), name=filt))

        # For a very small fraction of stars, the star phase falls on integers in-between
        # the ones we have definition for, as a result of the interpolation. For these
        # stars, round phase down to nearest defined phase (e.g., if phase is 71,
        # then round it down to 5, rather than up to 101).
        # Note: this only becomes relevant when the cluster is > 10**6 M-sun, this
        # effect is so small
        # Convert nan_to_num to avoid errors on greater than, less than comparisons
        star_systems_phase_non_nan = np.nan_to_num(star_systems['phase'], nan=-99)
        bad = (star_systems_phase_non_nan > 5) & (star_systems_phase_non_nan < 101) & (star_systems_phase_non_nan != 9) & (star_systems_phase_non_nan != -99)
        star_systems['phase'][bad] = 5


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
            good = r_id_tmp > 0
            idx_rem_good = idx_rem[good]

            star_systems['mass_current'][idx_rem_good] = r_mass_tmp[good]
            star_systems['phase'][idx_rem_good] = r_id_tmp[good]

            # Give remnants a magnitude of nan, so they can be filtered out later when calculating flux.
            for filt in self.filt_names:
                star_systems[filt][idx_rem_good] = np.full(len(idx_rem_good), np.nan)

        return star_systems


    def _make_companions_table_new(self, star_systems, compMass):
        """Make companions table for resolved clusters with multiplicity.

        Parameters
        ----------
        star_systems : astropy.table.Table
            Table containing the properties of the primary stars.
        compMass : numpy.ma.MaskedArray
            Masked array containing the masses of the companions.

        Returns
        -------
        companions : astropy.table.Table
        """
        N_systems = len(star_systems)
        N_companions = np.sum(~compMass.mask, axis=1)
        N_comp_tot = np.sum(N_companions)
        star_systems.add_column(Column(N_companions, name='N_companions'))
        system_index = np.repeat(np.arange(N_systems), N_companions)
        companions = Table([system_index], names=['system_idx'])
        companions.add_column(np.zeros(N_comp_tot, dtype=float), name='mass')

        if isinstance(self.imf._multi_props, multiplicity.MultiplicityResolvedDK):
            companions.add_column(Column(self.imf._multi_props.log_semimajoraxis(star_systems['mass'][companions['system_idx']]), name='log_a'))
            companions.add_column(Column(self.imf._multi_props.random_e(self.rng.random(N_comp_tot)), name='e'))
            companions['i'], companions['Omega'], companions['omega'] = self.imf._multi_props.random_keplarian_parameters(
                self.rng.random(N_comp_tot),
                self.rng.random(N_comp_tot),
                self.rng.random(N_comp_tot)
            )

        companions['mass'] = compMass.compressed()
        for key in ['Teff', 'L', 'logg', 'mass_current']:
            companions[key] = self.iso_interps[key](companions['mass'])

        for key in ['isWR', 'phase']:
            companions[key] = np.round(self.iso_interps[key](companions['mass']))

        companions['metallicity'] = np.ones(N_comp_tot) * self.iso.metallicity

        # For a very small fraction of stars, the star phase falls on integers in-between
        # the ones we have definition for, as a result of the interpolation. For these
        # stars, round phase down to nearest defined phase (e.g., if phase is 71,
        # then round it down to 5, rather than up to 101).
        # Convert nan_to_num to avoid errors on greater than, less than comparisons
        companions_phase_non_nan = np.nan_to_num(companions['phase'], nan=-99)
        companions['phase'][
            (companions_phase_non_nan > 5) &
            (companions_phase_non_nan < 101) &
            (companions_phase_non_nan != 9) &
            (companions_phase_non_nan != -99)
        ] = 5

        # Update primary fluxes to include the flux of companions.
        for filt in self.filt_names:
            companions[filt] = self.iso_interps[filt](companions['mass'])
            primary_flux = 10**(-star_systems[filt] / 2.5)
            # Sum the flux of all companions in each system
            companions_flux = np.bincount(companions['system_idx'], weights=10**(-companions[filt] / 2.5), minlength=N_systems)
            combined_flux = np.nansum(np.vstack((primary_flux, companions_flux)), axis=0)
            combined_flux[combined_flux == 0] = np.nan
            star_systems[filt] = -2.5 * np.log10(combined_flux)

        #####
        # Make Remnants with flux = 0 in all bands.
        #####
        if self.ifmr:
            # Identify compact objects as those with Teff = 0 or with masses above the max iso mass
            highest_mass_iso = self.iso.points['mass'].max()
            remnant_idx = np.where(np.isnan(companions['Teff']) & (companions['mass'] > highest_mass_iso))[0]
            self.remnant_idx_new = remnant_idx
            # Calculate remnant mass and ID for compact objects; update remnant_id and remnant_mass arrays accordingly
            if 'metallicity_array' in inspect.getfullargspec(self.ifmr.generate_death_mass).args:
                remnant_mass, remnant_code = self.ifmr.generate_death_mass(mass_array=companions['mass'][remnant_idx], metallicity_array=companions['metallicity'][remnant_idx])
            else:
                remnant_mass, remnant_code = self.ifmr.generate_death_mass(mass_array=companions['mass'][remnant_idx])

            # Drop remnants where it is not relevant (e.g. not a compact object or outside mass range IFMR is defined for)
            remnant_valid = remnant_code > 0
            remnant_valid_idx = remnant_idx[remnant_valid]
            self.remnant_mass_new = remnant_mass
            self.remnant_valid_idx_new = remnant_valid_idx
            companions['mass_current'][remnant_valid_idx] = remnant_mass[remnant_valid]
            companions['phase'][remnant_valid_idx] = remnant_code[remnant_valid]
            # Give remnants a magnitude of nan, so they can be filtered out later when calculating flux.
            for filt in self.filt_names:
                companions[filt][remnant_valid_idx] = np.full(len(remnant_idx[remnant_valid]), np.nan)

        companions_teff_non_nan = np.nan_to_num(companions['Teff'], nan=-99)
        if self.verbose and sum(companions_teff_non_nan > 0) != N_comp_tot:
            print(f'Found {N_comp_tot - sum(companions_teff_non_nan > 0):d} companions out of stellar mass range')

        assert companions['mass'][companions_teff_non_nan > 0].min() > 0, "Companion mass is not positive"

        return star_systems, companions


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

            companions['e'] = self.imf._multi_props.random_e(self.rng.random(N_comp_tot))
            companions['i'], companions['Omega'], companions['omega'] = self.imf._multi_props.random_keplarian_parameters(
                self.rng.random(N_comp_tot),
                self.rng.random(N_comp_tot),
                self.rng.random(N_comp_tot)
            )


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

            # companions['mass'][cdx] = compMass[idx, cc-1]
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
                companions_phase_non_nan = np.nan_to_num(companions['phase'], nan=-99)
                bad = np.where( (companions_phase_non_nan > 5) &
                                (companions_phase_non_nan < 101) &
                                (companions_phase_non_nan != 9) &
                                (companions_phase_non_nan != -99))
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
            cdx_rem = np.where(np.isnan(companions['Teff']) &
                                (companions['mass'] > highest_mass_iso))[0]

            # Calculate remnant mass and ID for compact objects; update remnant_id and
            # remnant_mass arrays accordingly
            if 'metallicity_array' in inspect.getfullargspec(self.ifmr.generate_death_mass).args:
                r_mass_tmp, r_id_tmp = self.ifmr.generate_death_mass(mass_array=companions['mass'][cdx_rem],
                                                                     metallicity_array=companions['metallicity'][cdx_rem])
            else:
                r_mass_tmp, r_id_tmp = self.ifmr.generate_death_mass(mass_array=companions['mass'][cdx_rem])


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
        # assert companions['mass'][idx].min() > 0, "Companion mass must be positive"

        return star_systems, companions


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
            idx = star_systems_teff_non_nan > 0
        else:
            # Keep stars (with Teff) and any other compact objects (with phase info).
            idx = (star_systems_teff_non_nan > 0) | (star_systems_phase_non_nan >= 0)

        if self.verbose and sum(idx) != N_systems:
            print( 'Found {0:d} stars out of mass range'.format(N_systems - sum(idx)))

        star_systems = star_systems[idx]
        N_systems = len(star_systems)

        if self.imf.make_multiples:
            # Clean up companion stuff (which we haven't handled yet)
            compMass = compMass[idx]

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
        rand_red = self.rng.standard_normal(len(self.star_systems))

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
                 rebin=True, verbose=False):

        self.verbose = verbose
        if self.verbose:
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
        L_all = 10**evol['logL'] * c.L_sun # luminsoity in W
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
        
        if self.verbose:
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
        The metallicity of the isochrone, in [M/H].
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
        the savefile exists. You should recompute anytime you change
        the filter set (see filters below).

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
                          'ubv,R', 'ubv,I'],
                verbose=False):
        self.metallicity = metallicity

        # Make the iso_dir, if it doesn't already exist
        if not os.path.exists(iso_dir):
            os.makedirs(iso_dir)

        # Make and input/output file name for the stored isochrone photometry.
        # For solar metallicity case, allow for legacy isochrones (which didn't have
        # metallicity tag since they were all solar metallicity) to be read
        # properly
        if metallicity == 0.0:
            save_file_fmt = '{0}/iso_{1:.2f}_{2:4.2f}_{3:4s}_p00.fits'
            self.save_file = save_file_fmt.format(iso_dir, logAge, AKs, str(round(distance)).zfill(5))

            save_file_legacy = '{0}/iso_{1:.2f}_{2:4.2f}_{3:4s}.fits'
            self.save_file_legacy = save_file_legacy.format(iso_dir, logAge, AKs, str(round(distance)).zfill(5))
        else:
            # Set metallicity flag
            if metallicity < 0:
                metal_pre = 'm'
            else:
                metal_pre = 'p'
            metal_flag = int(abs(metallicity)*10)

            save_file_fmt = '{0}/iso_{1:.2f}_{2:4.2f}_{3:4s}_{4}{5:2s}.fits'
            self.save_file = save_file_fmt.format(iso_dir, logAge, AKs, str(round(distance)).zfill(5), metal_pre, str(metal_flag).zfill(2))
            self.save_file_legacy = save_file_fmt.format(iso_dir, logAge, AKs, str(round(distance)).zfill(5), metal_pre, str(metal_flag).zfill(2))
        
        # Expected filters
        self.filters = filters

        # Recalculate isochrone if save_file doesn't exist or recomp == True
        file_exists = self.check_save_file(evo_model, atm_func, red_law, verbose=verbose)

        if (not file_exists) | (recomp==True):
            self.recalc = True
            if verbose:
                print(f'Generating new isochrone of log(t)={logAge:.2f}, AKs={AKs:.2f}, d={distance} pc')
            # user_input = input(f"Isochrone file {self.save_file} does not exist or needs to be regenerated. Do you want to proceed? (yes/no): ").strip().lower()
            # if user_input != 'yes':
            #     print("Operation canceled by the user.")
            #     return
            super().__init__(logAge, AKs, distance,
                             metallicity=metallicity,
                             evo_model=evo_model, atm_func=atm_func,
                             wd_atm_func=wd_atm_func,
                             wave_range=wave_range,
                             red_law=red_law, mass_sampling=mass_sampling,
                             min_mass=min_mass, max_mass=max_mass, rebin=rebin, verbose=verbose)

            # Make photometry
            self.make_photometry(rebin=rebin, vega=vega)
            if self.verbose:
                print(f'Isochrone saved to {self.save_file}')
        else:
            self.recalc = False
            try:
                self.points = Table.read(self.save_file)
            except:
                self.points = Table.read(self.save_file_legacy)
            # print(f'Isochrone loaded from existing file: {self.save_file}')
            # Add some error checking.

        return

    def make_photometry(self, rebin=True, vega=vega):
        """
        Make synthetic photometry for the specified filters. This function
        udpates the self.points table to include new columns with the
        photometry.

        """
        if self.verbose:
            startTime = time.time()

        # print('Making photometry for isochrone: log(t) = %.2f  AKs = %.2f  dist = %d' % \
        #     (meta['LOGAGE'], meta['AKS'], meta['DISTANCE']))
        # print( 'Starting at: ', datetime.datetime.now(), '  Usually takes ~5 minutes')

        npoints = len(self.points)
        verbose_fmt = 'M = {0:7.3f} Msun  T = {1:5.0f} K  m_{2:s} = {3:4.2f}'

        # Loop through the filters, get filter info, make photometry for
        # all stars in this filter.
        for ii in self.filters:
            if self.verbose:
                prt_fmt = 'Starting filter: {0:s}   Elapsed time: {1:.2f} seconds'
                print( prt_fmt.format(ii, time.time() - startTime))

            filt = get_filter_info(ii, rebin=rebin, vega=vega)
            filt_name = get_filter_col_name(ii)

            # Make the column to hold magnitudes in this filter. Add to points table.
            col_name = 'm_' + filt_name
            mag_col = Column(np.zeros(npoints, dtype=float), name=col_name)
            self.points.add_column(mag_col)

            # Loop through each star in the isochrone and do the filter integration
            if self.verbose:
                print('Starting synthetic photometry')
            for ss in range(npoints):
                star = self.spec_list[ss]  # These are already extincted, observed spectra.
                star_mag = mag_in_filter(star, filt)

                self.points[col_name][ss] = star_mag

                if self.verbose and (ss % 100) == 0:
                    print( verbose_fmt.format(self.points['mass'][ss], self.points['Teff'][ss],
                                             filt_name, star_mag))

        if self.verbose:
            endTime = time.time()
            print( 'Time taken: {0:.2f} seconds'.format(endTime - startTime))

        if self.save_file != None:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                self.points.write(self.save_file, overwrite=True)

        return

    def check_save_file(self, evo_model, atm_func, red_law, verbose=False):
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
            else:
                if verbose:
                    print(f'Isochrone file {self.save_file} exists, but meta-data does not match.')
                    if tmp.meta['EVOMODEL'] != type(evo_model).__name__:
                        print(f'  EVOMODEL: {tmp.meta["EVOMODEL"]} != {type(evo_model).__name__}')
                    if tmp.meta['ATMFUNC'] != atm_func.__name__:
                        print(f'  ATMFUNC: {tmp.meta["ATMFUNC"]} != {atm_func.__name__}')
                    if tmp.meta['REDLAW'] != red_law.name:
                        print(f'  REDLAW: {tmp.meta["REDLAW"]} != {red_law.name}')
        else:
            if verbose:
                print(f'Isochrone file {self.save_file} or {self.save_file_legacy} does not exist, generating new one.')

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

            gravity = float( logg_all[ii] )
            L = float( L_all[ii].cgs / (units.erg / units.s)) # in erg/s
            T = float( T_all[ii] / units.K)               # in Kelvin
            R = float( R_all[ii].to('pc') / units.pc)              # in pc

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
                    AKs_act = self.rng.normal(loc=AKs, scale=dAKs)
                    # Apply dAKs_max if desired. Redo if diff > dAKs_max
                    if dAKs_max != None:
                        diff = abs(AKs_act - AKs)
                        while diff > dAKs_max:
                            print('While loop active')
                            AKs_act = self.rng.normal(loc=AKs, scale=dAKs)
                            diff = abs(AKs_act - AKs)
                elif dist == 'uniform':
                    low = AKs - dAKs
                    high = AKs + dAKs
                    AKs_act = self.rng.uniform(low=low, high=high)
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

    elif name.startswith('rubin'):
        filt = filters.get_rubin_filt(filterName)

    else:
        # Otherwise, look for the filter info in the cdbs/mtab and cdbs/comp files
        try:
            filt = ObsBandpass(name)
        except:
            raise Exception('Filter {0} not understood. Check spelling and make sure cdbs/mtab and cdbs/comp files are up to date'.format(name))

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
        elif 'roman' in tmp:
            filt_name = 'roman_{0}'.format(tmp[-1])
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
                 'naco_IB_2.00': 'naco,IB_2.00', 'naco_IB_2.03':'naco,IB_2.03', 'naco_IB_2.06':'naco,IB_2.06',
                 'naco_IB_2.24':'naco,IB_2.24', 'naco_IB_2.27':'naco,IB_2.27',
                 'naco_IB_2.30':'naco,IB_2.30', 'naco_IB_2.33':'naco,IB_2.33',
                 'naco_IB_2.36':'naco,IB_2.36',
                 'ukirt_J':'ukirt,J', 'ukirt_H':'ukirt,H', 'ukirt_K':'ukirt,K',
                 'ctio_osiris_H': 'ctio_osiris,H', 'ctio_osiris_K': 'ctio_osiris,K',
                 'ztf_g':'ztf,g', 'ztf_r':'ztf,r', 'ztf_i':'ztf,i',
                 'gaiaDR2_G': 'gaia,dr2_rev,G', 'gaiaDR2_Gbp':'gaia,dr2_rev,Gbp',
                 'gaiaDR2_Grp':'gaia,dr2_rev,Grp',
                 'hawki_J': 'hawki,J',
                 'hawki_H': 'hawki,H',
                 'hawki_Ks': 'hawki,Ks',
                 'roman_f062': 'roman,wfi,f062',
                 'roman_f087': 'roman,wfi,f087',
                 'roman_f106': 'roman,wfi,f106',
                 'roman_f129': 'roman,wfi,f129',
                 'roman_f158': 'roman,wfi,f158',
                 'roman_w146': 'roman,wfi,w146',
                 'roman_f213': 'roman,wfi,f213',
                 'roman_f184': 'roman,wfi,f184',
                 'rubin_g':'rubin,g',
                 'rubin_i':'rubin,i',
                 'rubin_r':'rubin,r',
                 'rubin_u':'rubin,u',
                 'rubin_z':'rubin,z',
                 'rubin_y':'rubin,y'}

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

