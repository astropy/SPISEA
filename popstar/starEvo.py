#
#
#
import math
import logging
from numpy import searchsorted, genfromtxt

log = logging.getLogger('starEvo')


class StellarEvolutionGeneva():
    def __init__(self):
        r"""
        Define intrinsic properties for Geneva stellar models.
        """
        # define metallicity parameters for Geneva models
        self.z_solar = 0.02
        self.z_list = [0.01, 0.02, 0.03]
        self.z_file_map = {0.01: 'z01/', 0.02: 'z02/', 0.03: 'z03/'}
        
        # populate list of model masses (in solar masses)
        self.mass_list = [(0.1 + i*0.005) for i in range(181)]
        
        # populate list of isochrone ages (log scale)
        age_list = [round(5.5 + 0.01*i, 2) for i in range(190)]
        age_list += [round(7.4 + 0.05*i, 2) for i in range(12)]
        age_list += [round(math.log10(1.e8*x), 2) for x in range(1, 10)]
        age_list += [round(math.log10(1.e9*x), 2) for x in range(1, 10)]
        self.age_list = age_list
        
        # specify location of model files
        self.model_dir = '../models/geneva_merged/'
        
        
    def massTrack(self, mass=0.5, metallicity=0.0):
        r"""
        Extract an individual mass track from the Geneva collection.
        
        """
        
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Geneva collection.
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity
        
        # check age and metallicity are within bounds
        if (math.log10(age) < 5.5) or (math.log10(age) > 9.78):
            logger.error('Requested age is out of bounds.')
            
        if (z_defined < 0.01) or (z_defined > 0.03):
            logger.error('Requested metallicity is out of bounds.')
        
        # convert age (in yrs) to log scale and find nearest value in grid
        age_idx = searchsorted(self.age_list, math.log10(age), side='right')
        iso_file = 'iso_' + str(self.age_list[age_idx]) + '.dat'
        
        # find closest metallicity value
        z_idx = searchsorted(self.z_list, z_defined, side='left')
        z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file
        
        # return isochrone data
        return genfromtxt(full_iso_file, comments='#')
