import math
import logging
from numpy import searchsorted, genfromtxt
import numpy as np
import os
import glob
import pdb
import warnings
from astropy.table import Table
from popstar.utils import objects

log = logging.getLogger('evolution')

# Fetch root directory of evolution models.
try:
    models_dir = os.environ['POPSTAR_MODELS']
    models_dir += '/evolution/'
except KeyError:
    warnings.warn("POPSTAR_MODELS is undefined; functionality "
                  "will be SEVERELY crippled.")
    models_dir = ''
    
class StellarEvolution(object):
    def __init__(self, model_dir, age_list, mass_list, z_list):
        self.model_dir = model_dir
        self.z_list = z_list
        self.mass_list = mass_list
        self.age_list = age_list
        
        return

    def massTrack(self, mass, metallicity):
        pass

    def isochrone(self, age=1.e8, metallicity=0.0):
        pass
    
class Geneva(StellarEvolution):
    def __init__(self):
        r"""
        Define intrinsic properties for Geneva stellar models.
        """
        # populate list of model masses (in solar masses)
        mass_list = [(0.1 + i*0.005) for i in range(181)]
        
        # define metallicity parameters for Geneva models
        z_list = [0.01, 0.02, 0.03]
        
        # populate list of isochrone ages (log scale)
        age_list = [round(5.5 + 0.01*i, 2) for i in range(190)]
        age_list += [round(7.4 + 0.05*i, 2) for i in range(12)]
        age_list += [round(math.log10(1.e8*x), 2) for x in range(1, 10)]
        age_list += [round(math.log10(1.e9*x), 2) for x in range(1, 10)]
        age_list = age_list
        
        # specify location of model files
        model_dir = '../models/geneva_merged/'

        StellarEvolution.__init__(model_dir, age_list, mass_list, z_list)

        self.z_solar = 0.02
        self.z_file_map = {0.01: 'z01/', 0.02: 'z02/', 0.03: 'z03/'}
        
        
        
    def massTrack(self, mass=0.5, metallicity=0.0):
        r"""
        Extract an individual mass track from the Geneva collection.
        
        """
        return
        
    
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
        iso_file = 'iso_' + str(self.age_list[age_idx]) + '.fits'
        
        # find closest metallicity value
        z_idx = searchsorted(self.z_list, z_defined, side='left')
        z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file
        
        # return isochrone data
        return genfromtxt(full_iso_file, comments='#')

    def format_isochrones(input_iso_dir, output_iso_dir):
        return

#---------------------------------------#
# Now for the Ekstrom+12 Geneva models
#---------------------------------------#

class Ekstrom(StellarEvolution):
    def __init__(self):
        r"""
        Define intrinsic properties for the Ekstrom+12 Geneva stellar models.
        """
        # populate list of model masses (in solar masses)
        #mass_list = [(0.1 + i*0.005) for i in range(181)]
        
        # define metallicity parameters for Geneva models
        #z_list = [0.01, 0.02, 0.03]
        
        # populate list of isochrone ages (log scale)
        #age_list = [round(5.5 + 0.01*i, 2) for i in range(190)]
        #age_list += [round(7.4 + 0.05*i, 2) for i in range(12)]
        #age_list += [round(math.log10(1.e8*x), 2) for x in range(1, 10)]
        #age_list += [round(math.log10(1.e9*x), 2) for x in range(1, 10)]
        #age_list = age_list
        
        # specify location of model files
        #model_dir = '../models/geneva_merged/'

        #super().__init__(model_dir, age_list, mass_list, z_list)

        #self.z_solar = 0.02
        #self.z_file_map = {0.01: 'z01/', 0.02: 'z02/', 0.03: 'z03/'}
        
        
        
    def massTrack(self, mass=0.5, metallicity=0.0):
        r"""
        Extract an individual mass track from the Ekstrom+12 Geneva collection.
        
        """
        return
        
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Ekstrom+12 Geneva collection.
        """
        # convert metallicity to mass fraction
        #z_defined = self.z_solar*10.**metallicity
        
        # check age and metallicity are within bounds
        #if (math.log10(age) < 5.5) or (math.log10(age) > 9.78):
        #    logger.error('Requested age is out of bounds.')
            
        #if (z_defined < 0.01) or (z_defined > 0.03):
        #    logger.error('Requested metallicity is out of bounds.')
        
        # convert age (in yrs) to log scale and find nearest value in grid
        #age_idx = searchsorted(self.age_list, math.log10(age), side='right')
        #iso_file = 'iso_' + str(self.age_list[age_idx]) + '.fits'
        
        # find closest metallicity value
        #z_idx = searchsorted(self.z_list, z_defined, side='left')
        #z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        #full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file
        
        # return isochrone data
        return genfromtxt(full_iso_file, comments='#')

    def format_isochrones(input_iso_dir):
        r"""
        Parse iso.fits (filename hardcoded) file downloaded from Ekstrom+12
        models, create individual isochrone files for the different ages.

        input_iso_directory should lead to 
            Ekstrom2012/iso/<metallicity> 
        directory, where iso.fits file should be located.

        Creates two new directories, rot and norot, which contain their 
        respective isochrones.
        """
        # Store current directory for later
        start_dir = os.getcwd()

        # Move into metallicity direcotry, read iso.fits file
        os.chdir(input_iso_dir)
    
        print 'Read Input: this is slow'
        iso = Table.read('iso.fits')
        print 'Done'    
    
        ages_all = iso['col1']
    
        # Extract the unique ages
        age_arr = np.unique(ages_all)

        # For each unique age, extract the proper rows and make corresponding
        # table. Be sure to separate rotating from non-rotating, and put in
        # separate subdirectories

        # First make the rot and norot directories, if they don't exit
        if os.path.exists('rot'):
            pass
        else:
            os.mkdir('rot')
            os.mkdir('norot')
    
        print 'Making individual isochrone files'
        for age in age_arr:
            good = np.where(ages_all == age)

            # Identify rot vs. non-rot
            idx_r = np.where(iso[good]['col2'] == 'r')
            idx_n = np.where(iso[good]['col2'] == 'n')

            tmp_r = iso[good][idx_r]
            tmp_n = iso[good][idx_n]
        
            # Write tables
            tmp_r.write('rot/iso_{0:4.2f}.fits'.format(age))
            tmp_n.write('norot/iso_{0:4.2f}.fits'.format(age))

        # Return to starting directory
        os.chdir(start_dir)

        return

    def create_iso(fileList, ageList, rot=True):
        """
        Given a set of isochrone files downloaded from
        http://obswww.unige.ch/Recherche/evoldb/index/Isochrone/, put in correct
        iso.fits format for parse_iso code.

        fileList: list of downloaded isochrone files (could be one)
    
        ageList: list of lists of ages associated with each file in filelist.
        MUST BE IN SAME ORDER AS ISOCHRONES IN FILE! Also needs to be in logAge
    
        rot = TRUE: assumes that models are rotating, will add appropriate column
    
        This code writes the individual files, which is then easiest to combine by hand
        in aquamacs 
        """
        # Read each file in fileList individually, add necessary columns
        for i in range(len(fileList)):
            t = Table.read(fileList[i],format='ascii')
            ages = ageList[i]

            # Find places where new models start; mass here is assumed to be 0.8
            start = np.where(t['M_ini'] == 0.8)

            # Now, each identified start is assumed to be associated with the
            # corresponding age in ages        
            if len(start[0]) != len(ages):
                print 'Ages mismatched in file! Quitting...'
                return

            age_arr = np.zeros(len(t))

        
            for j in range(len(start[0])):
                low_ind = start[0][j]
                # Deal with case at end of file
                if (j == len(start[0])-1):
                    high_ind = len(t)
                else:
                    high_ind = start[0][j+1]

                ind = np.arange(low_ind, high_ind, 1)
                age_arr[ind] = ages[j]

            # Add ages_arr column to column 1 in ischrone, as well as column
            # signifying rotation
            col_age = Column(age_arr, name = 'logAge')
            rot_val = np.chararray(len(t))
            rot_val[:] = 'r'
            if not rot:
                rot_val[:] = 'n'
            
            col_rot = Column(rot_val, name='Rot')
        
            t.add_column(col_rot, index=0)
            t.add_column(col_age, index=0)

            t.write('tmp'+str(i)+'.fits')

        return

#---------------------------------------#
# Now for the Parsec version 1.2s models
#---------------------------------------#

class Parsec(StellarEvolution):
    def __init__(self):
        r"""
        Define intrinsic properties for the Parsec version 1.2s stellar
        models.
        """
        # populate list of model masses (in solar masses)
        #mass_list = [(0.1 + i*0.005) for i in range(181)]
        
        # define metallicity parameters for Geneva models
        #z_list = [0.01, 0.02, 0.03]
        
        # populate list of isochrone ages (log scale)
        #age_list = [round(5.5 + 0.01*i, 2) for i in range(190)]
        #age_list += [round(7.4 + 0.05*i, 2) for i in range(12)]
        #age_list += [round(math.log10(1.e8*x), 2) for x in range(1, 10)]
        #age_list += [round(math.log10(1.e9*x), 2) for x in range(1, 10)]
        #age_list = age_list
        
        # specify location of model files
        #model_dir = '../models/geneva_merged/'

        #super().__init__(model_dir, age_list, mass_list, z_list)

        #self.z_solar = 0.02
        #self.z_file_map = {0.01: 'z01/', 0.02: 'z02/', 0.03: 'z03/'}
        
        
        
    def massTrack(self, mass=0.5, metallicity=0.0):
        r"""
        Extract an individual mass track from the Parsec version 1.2s
        collection.
        
        """
        return
        
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Parsec version 1.2s
        collection.
        """
        # convert metallicity to mass fraction
        #z_defined = self.z_solar*10.**metallicity
        
        # check age and metallicity are within bounds
        #if (math.log10(age) < 5.5) or (math.log10(age) > 9.78):
        #    logger.error('Requested age is out of bounds.')
            
        #if (z_defined < 0.01) or (z_defined > 0.03):
        #    logger.error('Requested metallicity is out of bounds.')
        
        # convert age (in yrs) to log scale and find nearest value in grid
        #age_idx = searchsorted(self.age_list, math.log10(age), side='right')
        #iso_file = 'iso_' + str(self.age_list[age_idx]) + '.fits'
        
        # find closest metallicity value
        #z_idx = searchsorted(self.z_list, z_defined, side='left')
        #z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        #full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file
        
        # return isochrone data
        return genfromtxt(full_iso_file, comments='#')

    def format_isochrones(input_iso_dir, metallicity_list):
        r"""
        Parse isochrone file downloaded from Parsec version 1.2 for different
        metallicities, create individual isochrone files for the different ages.
    
        input_iso_dir: points to ParsecV1.2s/iso directory. Assumes metallicity
        subdirectories already exist with isochrone files downloaded in them
        (isochrones files expected to start with "output*")

        metallicity_list format: absolute (vs. relative to solar),
        z + <digits after decimal>: e.g. Z = 0.014 --> z014
        """
        # Store current directory for later
        start_dir = os.getcwd()

        # Move into isochrone directory
        os.chdir(input_iso_dir)
        
        # Work on each metallicity isochrones individually
        for metal in metallicity_list:
            # More into metallicity directory, read isochrone file
            os.chdir(metal)

            isoFile = glob.glob('output*')
            print 'Read Input: this is slow'
            iso = Table.read(isoFile[0], format='fits')
            print 'Done'
    
            ages_all = iso['col2']

            # Extract the unique ages
            age_arr = np.unique(ages_all)

            # For each unique age, extract the proper rows and make corresponding
            # table
            print 'Making individual isochrone files'
            for age in age_arr:
                good = np.where(ages_all == age)
                tmp = iso[good]

                #Write table
                tmp.write('iso_{0:4.2f}.fits'.format(age))

            # Move back into iso directory
            os.chdir('..')

        # Return to starting directory
        os.chdir(start_dir)
        return


#---------------------------------------#
# Now for the Pisa (Tognelli+11) models
#---------------------------------------#

class Pisa(StellarEvolution):
    def __init__(self):
        r"""
        Define intrinsic properties for the Pisa (Tognelli+11) stellar
        models.
        """
        # populate list of model masses (in solar masses)
        #mass_list = [(0.1 + i*0.005) for i in range(181)]
        
        # define metallicity parameters for Geneva models
        #z_list = [0.01, 0.02, 0.03]
        
        # populate list of isochrone ages (log scale)
        #age_list = [round(5.5 + 0.01*i, 2) for i in range(190)]
        #age_list += [round(7.4 + 0.05*i, 2) for i in range(12)]
        #age_list += [round(math.log10(1.e8*x), 2) for x in range(1, 10)]
        #age_list += [round(math.log10(1.e9*x), 2) for x in range(1, 10)]
        #age_list = age_list
        
        # specify location of model files
        #model_dir = '../models/geneva_merged/'

        #super().__init__(model_dir, age_list, mass_list, z_list)

        #self.z_solar = 0.02
        #self.z_file_map = {0.01: 'z01/', 0.02: 'z02/', 0.03: 'z03/'}
        
        
        
    def massTrack(self, mass=0.5, metallicity=0.0):
        r"""
        Extract an individual mass track from the Pisa (Tognelli+11)
        collection.
        
        """
        return
        
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Pisa (Tognelli+11)
        collection.
        """
        # convert metallicity to mass fraction
        #z_defined = self.z_solar*10.**metallicity
        
        # check age and metallicity are within bounds
        #if (math.log10(age) < 5.5) or (math.log10(age) > 9.78):
        #    logger.error('Requested age is out of bounds.')
            
        #if (z_defined < 0.01) or (z_defined > 0.03):
        #    logger.error('Requested metallicity is out of bounds.')
        
        # convert age (in yrs) to log scale and find nearest value in grid
        #age_idx = searchsorted(self.age_list, math.log10(age), side='right')
        #iso_file = 'iso_' + str(self.age_list[age_idx]) + '.fits'
        
        # find closest metallicity value
        #z_idx = searchsorted(self.z_list, z_defined, side='left')
        #z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        #full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file
        
        # return isochrone data
        return genfromtxt(full_iso_file, comments='#')

    def format_isochrones(input_iso_dir, metallicity_list):
        r"""
        Rename the isochrone files extracted from Pisa (Tognelli+11) to fit
        naming/directory scheme

        input_iso_dir: points to Pisa2011/iso directory. Individual
        metallicity directories with the downloaded isochrones are
        expected to already exist there
        
        metallicity_list is the list of metallicities on which function
        is to be run.

        format for metallicity_list : absolute (vs. relative to sun)
        'z' + <digits after decimal>, e.g Z = 0.015 --> z015.
        """
        # Store current directory for later
        start_dir = os.getcwd()

        # Move into isochrone directory
        os.chdir(input_iso_dir)
        # Work on each metallicity directory individually
        for metal in metallicity_list:
            # Move into directory, check to see if files are already formatted
            os.chdir(metal)

            if os.path.exists('iso_6.00.fits'):
                print 'Files in {0:s} already formatted'.format(metal)
            else:
                # Create a ReadMe with the original file names to preserve the
                # model details
        
                cmd = "ls *.FITS > ReadMe"
                os.system(cmd)

                # Collect all filenames in a list, rename files one
                # by one
                isoFile_list = glob.glob('*.FITS')
                for File in isoFile_list:
                    name = File.split('_')
                    # Extract iso age from filename
                    age = float(name[1][1:])
                    logAge = np.log10(age * 10**6)

                    cmd = "mv {0:s} iso_{1:4.2f}.fits".format(File, logAge)
                    os.system(cmd)

            # Return to overhead directory
            os.chdir('..')

        # Return to starting directory
        os.chdir(start_dir)
        return

    def make_isochrone_grid(metallicity=0.015):
        """
        Create isochrone grid of given metallicity with time sampling = 0.01
        in logAge (hardcoded). This interpolates the downloaded isochrones
        when necessary. Builds upon the online iscohrone grid.

        Note: format of metallicity is important. After decimal point, must match
        the format of the metallcity directory (i.e., 0.015 matches directory z015,
        while 0.0150 would not)
        """
        logAge_arr = np.arange(6.0, 8.0+0.005, 0.01)
    
        count = 0
        for logAge in logAge_arr:
            # Could interpolate using evolutionary tracks, but less accurate.
            make_isochrone_pisa_interp(logAge, metallicity=metallicity)

            count += 1
        
            print 'Done {0} of {1} models'.format(count, (len(logAge_arr)))

        return
    
class MergedPisaEkstromParsec(StellarEvolution):
    def __init__(self):
        """
        Define intrinsic properties for merged Pisa-Ekstrom-Parsec
        stellar models.
        """
        # populate list of model masses (in solar masses)
        mass_list = [(0.1 + i*0.005) for i in range(181)]
        
        # define metallicity parameters for Geneva models
        z_list = [0.015]
        
        # populate list of isochrone ages (log scale)
        age_list = np.arange(6.0, 8.001, 0.01).tolist()
        
        # specify location of model files
        model_dir = models_dir + 'merged/pisa_ekstrom_parsec/'
        StellarEvolution.__init__(self, model_dir, age_list, mass_list, z_list)
        self.z_solar = 0.015
        self.z_file_map = {0.015: 'z015/'}

        
    def massTrack(self, mass=0.5, metallicity=0.0):
        r"""
        Extract an individual mass track from the Geneva collection.
        
        """
        return
        
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Geneva collection.
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)
        
        # check age and metallicity are within bounds
        if (log_age < self.age_list[0]) or (log_age > self.age_list[-1]):
            logger.error('Requested age is out of bounds.')
            
        if not z_defined in self.z_list:
            logger.error('Requested metallicity is out of bounds.')
        
        # convert age (in yrs) to log scale and find nearest value in grid
        age_idx = searchsorted(self.age_list, log_age, side='right')
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])
        
        # find closest metallicity value
        z_idx = searchsorted(self.z_list, z_defined, side='left')
        z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        full_iso_file = self.model_dir + z_dir + iso_file
        
        # return isochrone data
        iso = Table.read(full_iso_file, format='fits')
        iso.rename_column('col1', 'mass')
        iso.rename_column('col2', 'logT')
        iso.rename_column('col3', 'logL')
        iso.rename_column('col4', 'logg')
        iso.rename_column('col5', 'logT_WR')
        iso.rename_column('col6', 'model_ref')

        iso.meta['log_age'] = log_age
        iso.meta['metallicity'] = metallicity
        
        return iso

def make_isochrone_pisa_interp(log_age, metallicity=0.015, 
                         tracks=None, test=False):
    """
    Read in a set of isochrones and generate an isochrone at log_age
    that is well sampled at the full range of masses.

    Puts isochrones is Pisa2011/iso/<metal>/
    """
    # If logage > 8.0, quit immediately...grid doesn't go that high
    if log_age > 8.0:
        print 'Age too high for Pisa grid (max logAge = 8.0)'
        return

    # Directory with where the isochrones will go (both downloaded and interpolated)
    rootDir = models_dir + '/Pisa2011/iso/'
    metSuffix = 'z' + str(metallicity).split('.')[-1]
    rootDir += metSuffix + '/'

    # Can we find the isochrone directory?
    if not os.path.exists(rootDir):
        print 'Failed to find Pisa PMS isochrones for metallicity = ' + metSuffix
        return

    # Check to see if isochrone at given age already exists. If so, quit
    if os.path.exists(rootDir+'iso_{0:3.2f}.fits'.format(log_age)):
        print 'Isochrone at logAge = {0:3.2f} already exists'.format(log_age)
        return
    
    # Name/directory for interpolated isochrone
    isoFile = rootDir+'iso_%3.2f.fits' % log_age
    outSuffix = '_%.2f' % (log_age)

    print '*** Generating Pisa isochrone for log t = %3.2f and Z = %.3f' % \
        (log_age, metallicity)

    print time.asctime(), 'Getting original Pisa isochrones.'
    iso = get_orig_pisa_isochrones(metallicity=metallicity)

    # First thing is to find the isochrones immediately above and below desired
    # age
    iso_log_ages = iso.log_ages
    tmp = np.append(iso_log_ages, log_age)

    # Find desired age in ordered sequence; isolate model younger and older
    tmp.sort()
    good = np.where(tmp == log_age)
    young_model_logage = tmp[good[0]-1]
    old_model_logage = tmp[good[0]+1]
    
    # Isolate younger/older isochrones
    young_ind = np.where(iso.log_ages == young_model_logage)
    old_ind = np.where(iso.log_ages == old_model_logage)

    young_iso = iso.isochrones[young_ind[0]]
    old_iso = iso.isochrones[old_ind[0]]

    # Need both younger and older model on same temperature grid for time
    # interpolation. Will adopt mass grid of whichever model is closer in time
    if abs(young_model_logage - log_age) <= abs(old_model_logage - log_age):
        # Use young model mass grid
        young_iso, old_iso = interpolate_iso_tempgrid(young_iso, old_iso)
        
    else:
        # Use old model mass grid
        old_iso, young_iso = interpolate_iso_tempgrid(old_iso, young_iso)

    # Now, can interpolate in time over the two models. Do this star by star.
    # Work in linear time here!!
    numStars = len(young_iso.M)
    
    interp_iso = Isochrone(log_age)
    interp_iso.log_Teff = np.zeros(numStars, dtype=float)
    interp_iso.log_L = np.zeros(numStars, dtype=float)
    interp_iso.log_g = np.zeros(numStars, dtype=float)
    interp_iso.M = young_iso.M # Since mass grids should already be matched
    
    for i in range(numStars):
        # Do interpolations in linear space
        model_ages = [10**young_model_logage[0], 10**old_model_logage[0]]
        target_age = 10**log_age
        #model_ages = [young_model_logage[0], old_model_logage[0]]
        #target_age = log_age
        
        # Build interpolation functions
        Teff_arr = [10**young_iso.log_Teff[i], 10**old_iso.log_Teff[i]]
        logL_arr = [10**young_iso.log_L[i], 10**old_iso.log_L[i]]
        logg_arr = [10**young_iso.log_g[i], 10**old_iso.log_g[i]]
        
        f_log_Teff = interpolate.interp1d(model_ages, Teff_arr, kind='linear')
        f_log_L = interpolate.interp1d(model_ages, logL_arr, kind='linear')
        f_log_g = interpolate.interp1d(model_ages, logg_arr, kind='linear')

        interp_iso.log_Teff[i] = np.log10(f_log_Teff(target_age))
        interp_iso.log_L[i] = np.log10(f_log_L(target_age))
        interp_iso.log_g[i] = np.log10(f_log_g(target_age))

    # If indicated, plot new isochrone along with originals it was interpolated
    # from
    if test:
        py.figure(1)
        py.clf()
        py.plot(interp_iso.log_Teff, interp_iso.log_L, 'k-', label = 'Interp')
        py.plot(young_iso.log_Teff, young_iso.log_L, 'b-',
                label = 'log Age = {0:3.2f}'.format(young_model_logage[0]))
        py.plot(old_iso.log_Teff, old_iso.log_L, 'r-',
                label = 'log Age = {0:3.2f}'.format(old_model_logage[0]))
        rng = py.axis()
        py.xlim(rng[1], rng[0])
        py.xlabel('log Teff')
        py.ylabel('log L')
        py.legend()
        py.title('Pisa 2011 Isochrone at log t = %.2f' % log_age)
        py.savefig(rootDir + 'plots/interp_isochrone_at' + outSuffix + '.png')
    
    print time.asctime(), 'Finished.'

    # Write output to file, MUST BE IN SAME ORDER AS ORIG FILES
    _out = open(isoFile, 'w')
    
    _out.write('%10s  %10s  %10s  %10s\n' % 
               ('# log L', 'log Teff', 'Mass', 'log g'))
    _out.write('%10s  %10s  %10s  %10s\n' % 
               ('# (Lsun)', '(Kelvin)', '(Msun)', '(cgs)'))

    for ii in range(len(interp_iso.M)):
        _out.write('%10.4f  %10.4f  %10.4f  %10.4f\n' %
                   (interp_iso.log_L[ii], interp_iso.log_Teff[ii], interp_iso.M[ii],
                    interp_iso.log_g[ii]))

    _out.close()

    return

def get_orig_pisa_isochrones(metallicity=0.015):
    """
    Helper code to get the original pisa isochrones at given metallicity.
    These are downloaded online
    """
    pms_dir = models_dir + '/Pisa2011/iso/iso_orig/'
    metSuffix = 'z' + str(metallicity).split('.')[-1]
    pms_dir += metSuffix + '/'

    if not os.path.exists(pms_dir):
        print 'Failed to find Siess PMS isochrones for metallicity = ' + metSuffix
        return
    
    # Collect the isochrones
    files = glob.glob(pms_dir + '*.dat')
    count = len(files)

    data = objects.DataHolder()

    data.isochrones = []
    data.log_ages = []
    
    # Extract useful params from isochrones
    for ff in range(len(files)):
        d = Table.read(files[ff], format='ascii')

        # Extract logAge from filename
        log_age = float(files[ff].split('_')[2][:-4])

        # Create an isochrone object   
        iso = Isochrone(log_age)
        iso.M = d['col3']
        iso.log_Teff = d['col2']
        iso.log_L = d['col1']

        # If a log g column exist, extract it. Otherwise, calculate
        # log g from T and L and add column at end
        if len(d.keys()) == 3:
            
            # Calculate log g from T and L
            L_sun = 3.8 * 10**33 #cgs
            SB_sig = 5.67 * 10**-5 #cgs
            M_sun = 2. * 10**33 #cgs
            G_const = 6.67 * 10**-8 #cgs
        
            radius = np.sqrt( (10**d['col1'] * L_sun) /
                          (4 * np.pi * SB_sig *  (10**d['col2'])**4) )
            g = (G_const * d['col3'] * M_sun) / radius**2


            iso.log_g = np.log10(g.astype(np.float))
        else:
            iso.log_g = d['col4']
        
        data.isochrones.append(iso)
        data.log_ages.append(log_age)

        # If it doesn't already exist, add a column with logg vals. This will
        # be appended at the end
        if len(d.keys()) == 3:
            logg_col = Column(iso.log_g, name = 'col4')
            d.add_column(logg_col, index=3)
            d.write(files[ff],format='ascii')
    data.log_ages = np.array(data.log_ages)

    # Resort so that everything is in order of increasing age
    sdx = data.log_ages.argsort()
    data.masses = data.log_ages[sdx]
    data.isochrones = [data.isochrones[ss] for ss in sdx]

    return data

class Isochrone(object):
    def __init__(self, log_age):
        self.log_age = log_age

        
class MergedPisaEkstromParsec_norot(StellarEvolution):
    def __init__(self):
        """
        Define intrinsic properties for merged Pisa-Ekstrom-Parsec
        stellar models.
        """
        # populate list of model masses (in solar masses)
        mass_list = [(0.1 + i*0.005) for i in range(181)]
        
        # define metallicity parameters for Geneva models
        z_list = [0.015]
        
        # populate list of isochrone ages (log scale)
        age_list = np.arange(6.0, 8.001, 0.01).tolist()
        
        # specify location of model files
        model_dir = models_dir + 'merged/pisa_ekstrom_parsec/norot/'
        StellarEvolution.__init__(self, model_dir, age_list, mass_list, z_list)
        self.z_solar = 0.015
        self.z_file_map = {0.015: 'z015/'}

        
    def massTrack(self, mass=0.5, metallicity=0.0):
        r"""
        Extract an individual mass track from the Geneva collection.
        
        """
        return
        
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Geneva collection.
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)
        
        # check age and metallicity are within bounds
        if (log_age < self.age_list[0]) or (log_age > self.age_list[-1]):
            logger.error('Requested age is out of bounds.')
            
        if not z_defined in self.z_list:
            logger.error('Requested metallicity is out of bounds.')
        
        # convert age (in yrs) to log scale and find nearest value in grid
        age_idx = searchsorted(self.age_list, log_age, side='right')
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])
        
        # find closest metallicity value
        z_idx = searchsorted(self.z_list, z_defined, side='left')
        z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        full_iso_file = self.model_dir + z_dir + iso_file
        
        # return isochrone data
        iso = Table.read(full_iso_file, format='fits')
        iso.rename_column('col1', 'mass')
        iso.rename_column('col2', 'logT')
        iso.rename_column('col3', 'logL')
        iso.rename_column('col4', 'logg')
        iso.rename_column('col5', 'logT_WR')
        iso.rename_column('col6', 'model_ref')

        iso.meta['log_age'] = log_age
        iso.meta['metallicity'] = metallicity
        
        return iso
