import math
import logging
from numpy import searchsorted, genfromtxt
import numpy as np
import os
import glob
import pdb
import warnings
from astropy.table import Table, vstack, Column
from scipy import interpolate
import pylab as py
from spisea.utils import objects
from astropy import constants as cs
from astropy import units
import spisea.IsochroneMakerReformattedVersionNew
from spisea.IsochroneMakerReformattedVersionNew import extractor, reformatter
un = units

logger = logging.getLogger('evolution')

# Fetch root directory of evolution models.
try:
    models_dir = os.environ['SPISEA_MODELS']
    models_dir += '/evolution/'
except KeyError:
    warnings.warn("SPISEA_MODELS is undefined; functionality "
                  "will be SEVERELY crippled.")
    models_dir = ''
    
class StellarEvolution(object):
    """
    Base Stellar evolution class.
    Parameters
    ----------
    model_dir: path
        Directory path to evolution model files
    age_list: list
        List of ages
    mass_list: list
        List of masses
    z_list: list
        List of metallicities
    """
    def __init__(self, model_dir, age_list, mass_list, z_list):
        self.model_dir = model_dir
        self.z_list = z_list
        self.mass_list = mass_list
        self.age_list = age_list
        
        return
    
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
        model_dir = models_dir + 'geneva/'

        StellarEvolution.__init__(self, model_dir, age_list, mass_list, z_list)

        self.z_solar = 0.02
        self.z_file_map = {0.01: 'z01/', 0.02: 'z02/', 0.03: 'z03/'}
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Geneva collection.
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity
        
        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
            
        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            logger.error('Requested metallicity {0} is out of bounds.'.format(z_defined))
        
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


#---------------------------------------#
# Now for the Ekstrom+12 Geneva models
#---------------------------------------#

class Ekstrom12(StellarEvolution):
    """
    Evolution models from 
    `Ekstrom et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012A%26A...537A.146E/abstract>`_.
    Downloaded from `website <http://obswww.unige.ch/Recherche/evoldb/index/Isochrone/>`_.
    Parameters
    ----------
    rot: boolean, optional
        If true, then use rotating Ekstrom models. Default is true.
    """
    def __init__(self, rot=True):
        # define metallicity parameters for Ekstrom+12 models
        self.z_list = [0.014]
        
        # populate list of isochrone ages (log scale)
        self.age_list = np.arange(6.0, 8.0+0.005, 0.01)
        
        # Specify location of model files
        self.model_dir = models_dir+'Ekstrom2012/'

        # Specifying metallicity
        self.z_solar = 0.014
        self.z_file_map = {0.014: 'z014/'}

        # Specify rotation or not
        self.rot = rot
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Ekstrom+12 Geneva collection.
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)
        
        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
            
        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            logger.error('Requested metallicity {0} is out of bounds.'.format(z_defined))
        
        # Find nearest age in grid to input grid
        if log_age != self.age_list[0]:
            age_idx = searchsorted(self.age_list, log_age, side='right')
        else:
            age_idx = searchsorted(self.age_list, log_age, side='left')
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])
        
        # find closest metallicity value
        z_idx = searchsorted(self.z_list, z_defined, side='left')
        z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        if self.rot:  
            full_iso_file = self.model_dir + 'iso/' + z_dir + 'rot/' + iso_file
        else:
            full_iso_file = self.model_dir + 'iso/' + z_dir + 'norot/' + iso_file
        
        # Return isochrone data
        iso = Table.read(full_iso_file, format='fits')
        iso.rename_column('col4', 'Z')
        iso.rename_column('col1', 'logAge')
        iso.rename_column('col3', 'mass')
        iso.rename_column('col6', 'mass_current')
        iso.rename_column('col7', 'logL')
        iso.rename_column('col8', 'logT')
        iso.rename_column('col22', 'logg')
        iso.rename_column('col9', 'logT_WR')

        # Add isWR column
        isWR = Column([False] * len(iso), name='isWR')
        idx_WR = np.where(iso['logT'] != iso['logT_WR'])
        isWR[idx_WR] = True
        iso.add_column(isWR)

        # Add a phase column... everything is just a star.
        iso.add_column( Column(np.ones(len(iso)), name = 'phase'))

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso

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
    
        print( 'Read Input: this is slow')
        iso = Table.read('iso.fits')
        print( 'Done'    )
    
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
    
        print( 'Making individual isochrone files')
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
        the server, put in correct
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
                print( 'Ages mismatched in file! Quitting...')
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
    """
    Evolution models from 
    `Bressan et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012MNRAS.427..127B/abstract>`_,
    version 1.2s.
    Downloaded from `here <http://stev.oapd.inaf.it/cgi-bin/cmd>_`
    Notes
    -----
    Evolution model parameters used in download:
    * n_Reimers parameter (mass loss on RGB) = 0.2
    * photometric system: HST/WFC3 IR channel
    * bolometric corrections OBC from Girardi+08, based on ATLAS9 ODFNEW models
    * Carbon star bolometric corrections from Aringer+09
    * no dust
    * no extinction
    * Chabrier+01 mass function
    """
    def __init__(self):
        r"""
        Define intrinsic properties for the Parsec version 1.2s stellar
        models.
        """
        # populate list of model masses (in solar masses)
        #mass_list = [(0.1 + i*0.005) for i in range(181)]
        
        # define metallicity parameters for Parsec models
        self.z_list = [0.005, 0.015, 0.04]
        
        # populate list of isochrone ages (log scale)
        self.age_list = np.arange(6.6, 10.12+0.005, 0.01)
        self.age_list = np.append(6.40, self.age_list)
        
        # Specify location of model files
        self.model_dir = models_dir+'ParsecV1.2s/'

        # Specifying metallicity
        self.z_solar = 0.015
        self.z_file_map = {0.005: 'z005/', 0.015: 'z015/', 0.04: 'z04/'}
        
        
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Parsec version 1.2s
        collection.
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)
        
        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
            
        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            logger.error('Requested metallicity {0} is out of bounds.'.format(z_defined))
        
        # Find nearest age in grid to input grid
        if log_age != self.age_list[0]:
            age_idx = searchsorted(self.age_list, log_age, side='right')
        else:
            age_idx = searchsorted(self.age_list, log_age, side='left')
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])
        
        # find closest metallicity value
        z_idx = searchsorted(self.z_list, z_defined, side='left')
        z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file
        
        # return isochrone data
        iso = Table.read(full_iso_file, format='fits')
        iso.rename_column('col1', 'Z')
        iso.rename_column('col2', 'logAge')
        iso.rename_column('col3', 'mass')
        iso.rename_column('col4', 'mass_current')
        iso.rename_column('col5', 'logL')
        iso.rename_column('col6', 'logT')
        iso.rename_column('col7', 'logg')
        iso.rename_column('col15', 'phase')
        iso['logT_WR'] = iso['logT']

        # Parsec doesn't identify WR stars, so identify all as "False"
        isWR = Column([False] * len(iso), name='isWR')
        iso.add_column(isWR)
        
        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso
        

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
            print( 'Read Input: this is slow')
            iso = Table.read(isoFile[0], format='fits')
            print( 'Done')
    
            ages_all = iso['col2']

            # Extract the unique ages
            age_arr = np.unique(ages_all)

            # For each unique age, extract the proper rows and make corresponding
            # table
            print( 'Making individual isochrone files')
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
    """
    Evolution models from 
    `Tognelli et al. 2011 <https://ui.adsabs.harvard.edu/abs/2011A%26A...533A.109T/abstract>`_.
    
    Downloaded `online <http://astro.df.unipi.it/stellar-models/index.php?m=1>`_
    Notes
    ------
    Parameters used in download:
    * Y = middle value of 3 provided (changes for different metallicities)
    * mixing length = 1.68
    * Deuterium fraction: 2*10^-5 for Z = 0.015, 0.03; 4*10^-4 for 0.005
    """
    def __init__(self):
        r"""
        Define intrinsic properties for the Pisa (Tognelli+11) stellar
        models.
        """
        # define metallicity parameters for Pisa models
        self.z_list = [0.015]
        
        # populate list of isochrone ages (log scale)
        self.age_list = np.arange(6.0, 8.01+0.005, 0.01)
        
        # Specify location of model files
        self.model_dir = models_dir+'Pisa2011/'

        # Specifying metallicity
        self.z_solar = 0.015
        self.z_file_map = {0.015: 'z015/'}
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Pisa (Tognelli+11)
        collection.
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)
        
        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
            
        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            logger.error('Requested metallicity {0} is out of bounds for evolution model. Available z-vals: {1}.'.format(z_defined, self.z_list))
        
        # Find nearest age in grid to input grid
        if log_age != self.age_list[0]:
            age_idx = searchsorted(self.age_list, log_age, side='right')
        else:
            age_idx = searchsorted(self.age_list, log_age, side='left')
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])
        
        # find closest metallicity value
        z_idx = searchsorted(self.z_list, z_defined, side='left')
        z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file
        
        # return isochrone data
        iso = Table.read(full_iso_file, format='fits')
        iso.rename_column('col1', 'logL')
        iso.rename_column('col2', 'logT')
        iso.rename_column('col3', 'mass')
        iso.rename_column('col4', 'logg')
        iso['logT_WR'] = iso['logT']

        # Pisa models are too low for WR phase, add WR column with all False
        isWR = Column([False] * len(iso), name='isWR')
        iso.add_column(isWR)

        # Add columns for current mass and phase. 
        iso.add_column( Column(np.zeros(len(iso)), name = 'phase'))
        iso.add_column( Column(iso['mass'], name = 'mass_current'))

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso

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
                print( 'Files in {0:s} already formatted'.format(metal))
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
        
            print( 'Done {0} of {1} models'.format(count, (len(logAge_arr))))

        return

#==============================#
# Baraffe+15 models
#==============================#
class Baraffe15(StellarEvolution):
    """
    Evolution models published in 
    `Baraffe et al. 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...577A..42B/abstract>`_.
    Downloaded from `BHAC15 site <http://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/BHAC15_tracks>`_.
    """
    def __init__(self):
        # define metallicity parameters for Baraffe models
        self.z_list = [0.015]
        
        # populate list of isochrone ages (log scale)
        self.age_list = np.arange(6.0, 8.0+0.005, 0.01)
        
        # Specify location of model files
        self.model_dir = models_dir+'Baraffe15/'

        # Specifying metallicity
        self.z_solar = 0.015
        self.z_file_map = {0.015: 'z015/'}
        
    
    def isochrone(self, age=5.e7, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Baraffe+15
        collection.
        """
       # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)
        
        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
            
        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            logger.error('Requested metallicity {0} is out of bounds.'.format(z_defined))
        
        # Find nearest age in grid to input grid
        if log_age != self.age_list[0]:
            age_idx = searchsorted(self.age_list, log_age, side='right')
        else:
            age_idx = searchsorted(self.age_list, log_age, side='left')
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])
        
        # find closest metallicity value
        z_idx = searchsorted(self.z_list, z_defined, side='left')
        z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file
        
        # Read isochrone, get in proper format
        iso = Table.read(full_iso_file, format='fits')
        iso.rename_column('Mass', 'mass')
        iso.rename_column('logG', 'logg')
        iso['logT'] = np.log10(iso['Teff'])
        
        # Pisa models are too low for WR phase, add WR column with all False
        iso['logT_WR'] = iso['logT']
        isWR = Column([False] * len(iso), name='isWR')
        iso.add_column(isWR)

        # Add columns for current mass and phase. 
        iso.add_column( Column(np.zeros(len(iso)), name = 'phase'))
        iso.add_column( Column(iso['mass'], name = 'mass_current'))

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso

    def tracks_to_isochrones(self, tracksFile):
        r"""
        Create isochrones at desired age sampling (6.0 < logAge < 8.0,
        steps of 0.01; hardcoded) from the Baraffe+15 tracks downloaded
        online. 
        tracksFile: tracks.dat file downloaded from Baraffe+15, with format
        modified to be read in python
        
        Writes isochrones in iso/ subdirectory off of work directory. Will
        create this subdirectory if it doesn't already exist
        """
        tracks = Table.read(tracksFile, format='ascii')

        age_arr = np.arange(6.0, 8.0+0.005, 0.01)
        #age_arr = [6.28]
        
        # Loop through the masses, interpolating track over time at each.
        # Resample track properties at hardcoded ages
        masses = np.unique(tracks['col1'])

        mass_interp = []
        age_interp = []
        Teff_interp = []
        logL_interp = []
        logG_interp = []
        print( 'Begin looping over masses')
        cnt=0
        for mass in masses:
            idx = np.where(tracks['col1'] == mass)
            tmp = tracks[idx]

            # First, extract Teff, logL, and logG, eliminating
            # duplicated inputs (these crash the interpolator)
            good_Teff = np.where( np.diff(tmp['col3']) != 0 )
            good_logG = np.where( np.diff(tmp['col5']) != 0 )
            good_logL = np.where( np.diff(tmp['col4']) != 0 )

            # Interpolate Teff, logL, and logG using linear interpolator
            tck_Teff = interpolate.interp1d(tmp['col2'], tmp['col3'])
            tck_logL = interpolate.interp1d(tmp['col2'], tmp['col4'])
            tck_logG = interpolate.interp1d(tmp['col2'], tmp['col5'])            


            Teff = tck_Teff(age_arr)
            logL = tck_logL(age_arr)
            logG = tck_logG(age_arr)
            
            # Test interpolation if desired
            test=False
            if test:
                py.figure(1, figsize=(10,10))
                py.clf()
                py.plot(tmp['col2'], tmp['col3'], 'k.', ms=8)
                py.plot(age_arr, Teff, 'r-', linewidth=2)
                py.xlabel('logAge')
                py.ylabel('Teff')
                py.savefig('test_Teff.png')

                py.figure(2, figsize=(10,10))
                py.clf()
                py.plot(tmp['col2'], tmp['col4'], 'k.', ms=8)
                py.plot(age_arr, logL, 'r-', linewidth=2)
                py.xlabel('logAge')
                py.ylabel('logL')
                py.savefig('test_logL.png')

                py.figure(3, figsize=(10,10))
                py.clf()
                py.plot(tmp['col2'], tmp['col5'], 'k.', ms=8)
                py.plot(age_arr, logG, 'r-', linewidth=2)
                py.xlabel('logAge')
                py.ylabel('logG')
                py.savefig('test_logG.png')
                
                pdb.set_trace()
           
            # Build upon arrays of interpolated values
            mass_interp = np.concatenate((mass_interp, np.ones(len(Teff)) * mass))
            age_interp = np.concatenate((age_interp, age_arr))
            Teff_interp = np.concatenate((Teff_interp, Teff))
            logL_interp = np.concatenate((logL_interp, logL))
            logG_interp = np.concatenate((logG_interp, logG))

            print( 'Done {0} of {1}'.format(cnt, len(masses)))
            cnt+=1

        # Now, construct the iso_*.fits files for each age, write files
        # to iso subdirectory
        # First check to see if subdirectory exists
        if not os.path.exists('iso/'):
            os.mkdir('iso')

        # Now for the loop
        ages = np.unique(age_interp)
        print( 'Writing iso files')
        for age in ages:
            good = np.where( age_interp == age)

            t = Table( (mass_interp[good], Teff_interp[good], logL_interp[good],
                       logG_interp[good]), names=('Mass', 'Teff', 'logL', 'logG') )

            # Write out as fits table
            name = 'iso_{0:3.2f}.fits'.format(age)
            t.write('iso/'+name, format='fits', overwrite=True)

        return

    def test_age_interp(self, onlineIso, interpIso):
        r"""
        Compare one of our interpolated ischrones with one
        of the isochrones provided online by Baraffe+15. 
        """
        true_iso = Table.read(onlineIso, format='ascii')
        our_iso = Table.read(interpIso, format='fits')
        
        # Compare the two isochrones using plots. Look at mass vs. Teff,
        # mass vs. logG, mass vs. logL. Ideally these isochrones should
        # be identical
        py.figure(1, figsize=(10,10))
        py.clf()
        py.plot(true_iso['col1'], true_iso['col2'], 'k.', ms = 10)
        py.plot(our_iso['Mass'], our_iso['Teff'], 'r.', ms = 10)
        py.xlabel('Mass')
        py.ylabel('Teff')
        py.savefig('interp_test1.png')

        py.figure(2, figsize=(10,10))
        py.clf()
        py.plot(true_iso['col1'], true_iso['col3'], 'k.', ms = 10)
        py.plot(our_iso['Mass'], our_iso['logL'], 'r.', ms = 10)
        py.xlabel('Mass')
        py.ylabel('logL')
        py.savefig('interp_test2.png')

        py.figure(3, figsize=(10,10))
        py.clf()
        py.plot(true_iso['col1'], true_iso['col4'], 'k.', ms = 10)
        py.plot(our_iso['Mass'], our_iso['logG'], 'r.', ms = 10)
        py.xlabel('Mass')
        py.ylabel('logG')
        py.savefig('interp_test3.png')

        # Look at the difference between values (assumes the masses are lined up)
        Teff_diff = np.mean(abs(true_iso['col2'][7:] - our_iso['Teff']))
        logL_diff = np.mean(abs(true_iso['col3'][7:] - our_iso['logL']))
        logG_diff = np.mean(abs(true_iso['col4'][7:] - our_iso['logG']))
    
        print( 'Average abs difference in Teff: {0}'.format(Teff_diff))
        print( 'Average abs difference in logL: {0}'.format(logL_diff))
        print( 'Average abs difference in logg: {0}'.format(logG_diff))

        return

def compare_Baraffe_Pisa(BaraffeIso, PisaIso):
    """
    Compare the Baraffe isochrones to the Pisa isochrones, since they overlap
    over a significant portion of mass space.
    """
    b = Table.read(BaraffeIso, format='fits')
    p = Table.read(PisaIso, format='ascii')

    name = BaraffeIso.split('_')
    age = name[1][:4]
    
    # Extract paramters we need
    b_mass = b['Mass']
    b_logT = np.log10(b['Teff'])
    b_logL = b['logL']
    b_logG = b['logG']

    p_mass = p['col3']
    p_logT = p['col2']
    p_logL = p['col1']
    p_logG = p['col4']

    m05_b = np.where( abs(b_mass - 0.5) == min(abs(b_mass - 0.5)) )
    m05_p = np.where( abs(p_mass - 0.5) == min(abs(p_mass - 0.5)) )
    
    # Comparison plots
    py.figure(1, figsize=(10,10))
    py.clf()
    py.plot(b_logT, b_logL, 'k-', linewidth=2, label='Baraffe+15')
    py.plot(b_logT[m05_b], b_logL[m05_b], 'k.', ms=10)
    py.plot(p_logT, p_logL, 'r', linewidth=2, label='Pisa')
    py.plot(p_logT[m05_p], p_logL[m05_p], 'r.', ms=10)
    py.xlabel('logT')
    py.ylabel('logL')
    py.title(age)
    py.axis([4.4, 3.4, -3, 4])
    #py.gca().invert_xaxis()
    py.legend()
    py.savefig('BaraffePisa_comp_{0}.png'.format(age))

    py.figure(2, figsize=(10,10))
    py.clf()
    py.plot(b_mass, b_logL, 'k-', linewidth=2, label='Baraffe+15')
    py.plot(b_mass[m05_b], b_logL[m05_b], 'k.', ms=10)
    py.plot(p_mass, p_logL, 'r', linewidth=2, label='Pisa')
    py.plot(p_mass[m05_p], p_logL[m05_p], 'r.', ms=10)
    py.xlabel('Mass')
    py.ylabel('logL')
    py.title(age)
    #py.axis([4.4, 3.4, -3, 4])
    #py.gca().invert_xaxis()
    py.legend()
    py.savefig('BaraffePisa_comp_mass_{0}.png'.format(age))    

    return

#===============================#
# MIST v.1 (Choi+16)
#===============================#
class MISTv1(StellarEvolution):
    """
    Define intrinsic properties for the MIST v1 stellar
    models. 
    Models originally downloaded from `online server <http://waps.cfa.harvard.edu/MIST/interp_isos.html>`_.
    Parameters
    ----------
    version: '1.0' or '1.2', optional
        Specify which version of MIST models you want. Version 1.0
        was downloaded from MIST website on 2/2017, while Version 1.2
        was downloaded on 8/2018 (solar metallicity)
        and 4/2019 (other metallicities). Default is 1.2.
    """
    def __init__(self, version=1.2):
        # define metallicity parameters for Parsec models
        self.z_list = [0.0000015,   # [Fe/H] = -4.00
                       0.0000047,   # [Fe/H] = -3.50
                       0.000015,    # [Fe/H] = -3.00
                       0.000047,    # [Fe/H] = -2.50
                       0.00015, # [Fe/H] = -2.00
                       0.00026, # [Fe/H] = -1.75
                       0.00047, # [Fe/H] = -1.50
                       0.00084, # [Fe/H] = -1.25
                       0.0015,  # [Fe/H] = -1.00
                       0.0026,  # [Fe/H] = -0.75
                       0.0046,  # [Fe/H] = -0.50
                       0.0082,  # [Fe/H] = -0.25
                       0.015,   # [Fe/H] = 0.00
                       0.024,   # [Fe/H] = 0.25
                       0.041]   # [Fe/H] = 0.50
        
        # populate list of isochrone ages (log scale)
        self.age_list = np.arange(5.01, 10.30+0.005, 0.01)

        # Set version directory
        self.version = version
        if self.version == 1.0:
            version_dir = 'v1.0/'
        elif self.version == 1.2:
            version_dir = 'v1.2/'
        else:
            raise ValueError('Version {0} not supported for MIST isochrones'.format(version))
        
        # Specify location of model files
        self.model_dir = models_dir+'MISTv1/' + version_dir

        # Specifying metallicity
        self.z_solar = 0.0142
        self.z_file_map = {0.0000015: 'z0000015/',
                           0.0000047: 'z0000047/',
                           0.000015: 'z000015/',
                           0.000047: 'z000047/',
                           0.00015: 'z00015/',
                           0.00026: 'z00026/',
                           0.00047: 'z00047/',
                           0.00084: 'z00084/',
                           0.0015: 'z0015/',
                           0.0026: 'z0026/',
                           0.0046: 'z0046/',
                           0.0082: 'z0082/',
                           0.015: 'z015/',
                           0.024: 'z024/',
                           0.041: 'z041/'}
        
        
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the MISTv1
        collection.
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar * (10.**metallicity)

        log_age = math.log10(age)

        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
            
        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            logger.error('Requested metallicity {0} is out of bounds.'.format(z_defined))

        # Find nearest age in grid to input grid
        if log_age != self.age_list[0]:
            age_idx = searchsorted(self.age_list, log_age, side='right')
        else:
            age_idx = searchsorted(self.age_list, log_age, side='left')
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])
        
        # find closest metallicity value
        z_idx = searchsorted(self.z_list, z_defined, side='left')
        if z_idx == len(self.z_list):   # in case just over last index
            z_idx = z_idx - 1
        z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file
        
        # return isochrone data. Column locations depend on
        # version
        iso = Table.read(full_iso_file, format='fits')
        if self.version == 1.0:
            iso.rename_column('col7', 'Z')
            iso.rename_column('col2', 'logAge')
            iso.rename_column('col3', 'mass')
            iso.rename_column('col4', 'logT')
            iso.rename_column('col5', 'logg')
            iso.rename_column('col6', 'logL')
            iso.rename_column('col65', 'phase')
        elif self.version == 1.2:
            iso.rename_column('col2', 'logAge')
            iso.rename_column('col3', 'mass')
            iso.rename_column('col4', 'mass_current')
            iso.rename_column('col9', 'logL')
            iso.rename_column('col14', 'logT')
            iso.rename_column('col17', 'logg')
            iso.rename_column('col79', 'phase')

        # For MIST isochrones, anything with phase = 6 is a WD.
        # Following our IFMR convention, change the phase designation
        # to 101
        isWD = np.where(iso['phase'] == 6)[0]
        iso['phase'][isWD] = 101

        # Define "isWR" column based on phase info
        isWR = Column([False] * len(iso), name='isWR')
        idx_WR = np.where(iso['phase'] == 9)[0]
        isWR[idx_WR] = True
        iso.add_column(isWR)

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso
        
    def format_isochrones(self, input_iso_dir, metallicity_list):
        r"""
        Parse isochrone file downloaded from MIST web server,
        create individual isochrone files for the different ages.
        Assumes all files start with MIST_iso*
        Parameters:
        -----------
        input_iso_dir: path
            Points to MISTv1/<version>/iso directory.
        metallicity_list: array
            List of metallicity directories to check (i.e. z015 is solar)
        """
        # Store current directory for later
        start_dir = os.getcwd()

        # Move into isochrone directory
        os.chdir(input_iso_dir)
        
        # Work on each metallicity isochrones individually
        for metal in metallicity_list:
            # More into metallicity directory, read isochrone file
            os.chdir(metal)

            isoFile = glob.glob('MIST_iso*')
            print( 'Read Input: this is slow')
            iso_f = Table()
            for ii in isoFile:
                tmp = Table.read(ii, format='ascii')
                iso_f = vstack([iso_f, tmp])
            print( 'Done')

            # Extract the unique ages
            ages_all = iso_f['col2']
            age_arr = np.unique(ages_all)

            # For each unique age, extract the proper rows and make corresponding
            # table
            print( 'Making individual isochrone files')
            for age in age_arr:
                good = np.where(ages_all == age)
                tmp = iso_f[good]

                # Need to make sure the tables are unmasked...this causes
                # problems later
                tmp2 = Table(tmp, masked=False)

                #Write table
                tmp2.write('iso_{0:4.2f}.fits'.format(age))

            # Move back into iso directory
            os.chdir('..')

        # Return to starting directory
        os.chdir(start_dir)
        return

#==============================#
# Merged model classes
#==============================#
class MergedBaraffePisaEkstromParsec(StellarEvolution):
    """
    This is a combination of several different evolution models:
    * Baraffe (`Baraffe et al. 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...577A..42B/abstract>`_)
    * Pisa (`Tognelli et al. 2011 <https://ui.adsabs.harvard.edu/abs/2011A%26A...533A.109T/abstract>`_)
    * Geneva (`Ekstrom et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012A%26A...537A.146E/abstract>`_)
    * Parsec (version 1.2s, `Bressan+12 <https://ui.adsabs.harvard.edu/abs/2012MNRAS.427..127B/abstract>`_)
    The model used depends on the age of the population and what stellar masses
    are being modeled:
    
    For logAge < 7.4:
    * Baraffe: 0.08 - 0.4 M_sun
    * Baraffe/Pisa transition: 0.4 - 0.5 M_sun 
    * Pisa: 0.5 M_sun to the highest mass in Pisa isochrone (typically 5 - 7 Msun)
    * Geneva: Highest mass of Pisa models to 120 M_sun
    For logAge > 7.4:
    * Parsec v1.2s: full mass range
    
    Parameters
    ----------
    rot: boolean, optional
        If true, then use rotating Ekstrom models. Default is true.
    """
    def __init__(self, rot=True):
        # populate list of model masses (in solar masses)
        mass_list = [(0.1 + i*0.005) for i in range(181)]
        
        # define metallicity parameters for Geneva models
        z_list = [0.015]
        
        # populate list of isochrone ages (log scale)
        age_list = np.arange(6.0, 10.091, 0.01).tolist()
        
        # specify location of model files
        model_dir = models_dir + 'merged/baraffe_pisa_ekstrom_parsec/'
        StellarEvolution.__init__(self, model_dir, age_list, mass_list, z_list)
        self.z_solar = 0.015
        
        # Switch to specify rotating/non-rotating models
        if rot:
            self.z_file_map = {0.015: 'z015_rot/'}
        else:
            self.z_file_map = {0.015: 'z015_norot/'}
        
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Baraffe-Pisa-Ekstrom-Parsec 
        collection
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)
        
        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
            
        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            logger.error('Requested metallicity {0} is out of bounds.'.format(z_defined))
        
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
        iso.rename_column('col6', 'mass_current')
        iso.rename_column('col7', 'phase')
        iso.rename_column('col8', 'model_ref')

        # Define "isWR" column based on phase info
        isWR = Column([False] * len(iso), name='isWR')
        idx_WR = np.where(iso['logT'] != iso['logT_WR'])
        isWR[idx_WR] = True
        iso.add_column(isWR)
        
        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)
        
        return iso


class MergedPisaEkstromParsec(StellarEvolution):
    """
    Same as MergedBaraffePisaEkstromParsec, but without
    the Baraffe models. 
    Parameters
    ----------
    rot: boolean, optional
        If true, then use rotating Ekstrom models. Default is true.
    """
    def __init__(self, rot=True):
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

        #Switch to specify rot/notot
        if rot:
            self.z_file_map = {0.015: 'z015_rot/'}
        else:
            self.z_file_map = {0.015: 'z015_norot/'}
        
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Pisa-Ekstrom-Parsec collection.
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)
        
        # check age and metallicity are within bounds
        if (log_age < self.age_list[0]) or (log_age > self.age_list[-1]):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
            
        if not z_defined in self.z_list:
            logger.error('Requested metallicity {0} is out of bounds.'.format(z_defined))
        
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
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)
        
        return iso

class MergedSiessGenevaPadova(StellarEvolution):
    """
    This is a combination of several different evolution models.
    The model used depends on the age of the population and what stellar masses
    are being modeled:
    * Siess (`Siess et al. 2000 <https://ui.adsabs.harvard.edu/abs/2000A%26A...358..593S/abstractt>`_)
    * Geneva (`Meynet & Maeder 2003 <https://ui.adsabs.harvard.edu/abs/2003A%26A...404..975M/abstract>`_)
    * Padova (`Marigo et al. 2008 <https://ui.adsabs.harvard.edu/abs/2008A%26A...482..883M/abstract>`_)
    For logAge < 7.4:
    
    * Siess: 0.1 - 7 M_sun
    * Siess/Geneva transition: 7 - 9 M_sun
    * Geneva: > 9 M_sun
    For logAge > 7.4:
    
    * Padova: full mass range
    """
    def __init__(self):
        """
        Define intrinsic properties for merged Siess-meynetMaeder-Padova 
        stellar models.
        """
        # populate list of model masses (in solar masses)
        mass_list = [(0.1 + i*0.005) for i in range(181)]
        
        # define metallicity parameters for Geneva models
        z_list = [0.02]
        
        # populate list of isochrone ages (log scale)
        age_list = np.arange(5.5, 7.41, 0.01).tolist()
        age_list.append(7.48)
        idx = np.arange(7.50, 8.01, 0.05)
        for ii in idx:
            age_list.append(ii)
        age_list.append(8.30)
        age_list.append(8.48)
        age_list.append(8.60)
        age_list.append(8.70)
        age_list.append(8.78)
        age_list.append(8.85)
        age_list.append(8.90)
        age_list.append(8.95)
        age_list.append(9.00)
        age_list.append(9.30)
        age_list.append(9.60)
        age_list.append(9.70)
        age_list.append(9.78)
        
        # specify location of model files
        model_dir = models_dir + 'merged/siess_meynetMaeder_padova/'
        StellarEvolution.__init__(self, model_dir, age_list, mass_list, z_list)
        self.z_solar = 0.02
        
        # Metallicity map
        self.z_file_map = {0.02: 'z02/'}
        
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Siess-Geneva-Padova collection.
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)
        
        # check age and metallicity are within bounds
        if (log_age < self.age_list[0]) or (log_age > self.age_list[-1]):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
            
        if not z_defined in self.z_list:
            logger.error('Requested metallicity {0} is out of bounds.'.format(z_defined))
        
        # convert age (in yrs) to log scale and find nearest value in grid
        age_idx = searchsorted(self.age_list, log_age, side='right')
        iso_file = 'iso_{0:.2f}.dat'.format(self.age_list[age_idx])
        
        # find closest metallicity value
        z_idx = searchsorted(self.z_list, z_defined, side='left')
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # generate isochrone file string
        full_iso_file = self.model_dir + z_dir + iso_file

        # return isochrone data
        iso = Table.read(full_iso_file, format='ascii')
        iso.rename_column('col1', 'mass')
        iso.rename_column('col2', 'logT')
        iso.rename_column('col3', 'logL')
        iso.rename_column('col4', 'logg')
        iso.rename_column('col5', 'logT_WR')
        iso.rename_column('col6', 'model_ref')
        
        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)
        
        return iso

#================================================#
    
def make_isochrone_pisa_interp(log_age, metallicity=0.015, 
                         tracks=None, test=False):
    """
    Read in a set of isochrones and generate an isochrone at log_age
    that is well sampled at the full range of masses.
    Puts isochrones is Pisa2011/iso/<metal>/
    """
    # If logage > 8.0, quit immediately...grid doesn't go that high
    if log_age > 8.0:
        print( 'Age too high for Pisa grid (max logAge = 8.0)')
        return

    # Directory with where the isochrones will go (both downloaded and interpolated)
    rootDir = models_dir + '/Pisa2011/iso/'
    metSuffix = 'z' + str(metallicity).split('.')[-1]
    rootDir += metSuffix + '/'

    # Can we find the isochrone directory?
    if not os.path.exists(rootDir):
        print( 'Failed to find Pisa PMS isochrones for metallicity = ' + metSuffix)
        return

    # Check to see if isochrone at given age already exists. If so, quit
    if os.path.exists(rootDir+'iso_{0:3.2f}.fits'.format(log_age)):
        print( 'Isochrone at logAge = {0:3.2f} already exists'.format(log_age))
        return
    
    # Name/directory for interpolated isochrone
    isoFile = rootDir+'iso_%3.2f.fits' % log_age
    outSuffix = '_%.2f' % (log_age)

    print( '*** Generating Pisa isochrone for log t = %3.2f and Z = %.3f' % \
        (log_age, metallicity))

    print( time.asctime(), 'Getting original Pisa isochrones.')
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
    
    print( time.asctime(), 'Finished.')

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
        print( 'Failed to find Siess PMS isochrones for metallicity = ' + metSuffix)
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

# ============================== #
# BPASS v2.2 models
# ============================== #
# These models were published by Stanway and Elridge et. al in 2018
# The code is based from the Merged Evolution model methods outlined in the
# existing code.
# After all, I want some degree of consistency with regards to
# Our original isochrone files are the stellar models (single, binary, and
# secondary) listed in the BPASS input files from one of the stellar models
# for ages 10**6.0, 10**6.1, ... , 10**11.0 years.


class BPASS(StellarEvolution):
    """ These models were published by Stanway and Elridge et. al in 2018
    Some of the code is based from the Merged Evolution model methods
    by Hosek et. al. After all, I want some degree of consistency with
    regards to how isochrones are made.Also, I do not want to reinvent
    the wheel twice.
    If one has not checked the IsochroneMakerReformatted.py file, each
    isochrone file is made of star from stellar model that are closest
    to specified input age and are within the specified margin of error
    for log(Age in years). This is to be doubly sure that I am getting
    the same age.
    Brief description: Creates a BPASS evolution object and
    specifies valid primary masses,
    secondary masses,  metallicities, ages.
    """
    def __init__(self):
        """
        Note: All of the possible mass list items include possible mass for
        single star, binary primary and secondary.
        Recall that BPASS age bins
        are 10**6.0, 10**6.1, ... , 10**11.0 yrs """
        self.age_list = [round(6.0 + x * 0.1, 1) for x in range(0, 51)]
        self.z_list = [
            10 ** -5,
            10 ** -4,
            10 ** -3,
            0.002,
            0.003,
            0.004,
            0.006,
            0.008,
            0.010,
            0.014,
            0.020,
            0.030,
            0.040,
            ]
        self.models_dir = models_dir + "/BPASS/v2.2/"
        self.z_solar = 0.020
        # The possible Metallicities of a BPASS Stellar model
        self.z_file_map = {
            10 ** -5: 'zem5/',
            10 ** -4: 'zem4/',
            0.00300: 'z003/',
            0.00200: 'z002/',
            0.00300: 'z003/',
            0.00400: 'z004/',
            0.00600: 'z006/',
            0.00800: 'z008/',
            0.01000: 'z010/',
            0.01400: 'z014/',
            0.02000: 'z020/',
            0.03000: 'z030/',
            0.04000: 'z040/',
            }
    def run_formatter(self, destination="/g/lu/scratch/ryotainagaki/" +
                      "BPASS_tester_newReformatTest/",
                      metallicity=["zem5", "zem4", "z001", "z002", "z003", "z004",
                                   "z006", "z008", "z010", "z014", "z030", "z040"]):
        """
        Reformat all BPASS stellar system models in  HOKI.MODELS_PATH
        using the reformatter function into fits tables with the following
        naming convention:
        Fitsmodels<original name of BPASS stellar system model file>.fits
        The files will be stored in destination/<BPASS metallicity> directory
        """
        reformatter(destination, metallicity)
        return
    
    def run_extractor(self, source="/g/lu/scratch/ryotainagaki/" +
                      "BPASS_tester_newReformatTest/",
                      metallicity=["zem5", "zem4", "z001", "z002", "z003",
                                   "z004", "z006", "z008", "z010", "z014",
                                   "z030", "z040"],
                     destination=models_dir,
                     times=51):
        """
        Uses the reformatted BPASS stellar system model files in source directory
        to create isochrones for log10(Age) = 6.0, 6.1, ... , 6.0 + 0.1*(x-times)
        and for every metallicity. Stores those isochrones into the
        destination/iso/<metallicity> directory.
        """
        for met in mets:
            for x in range(times):
                extractor(round(6.0 + 0.1 * x,1), met, source,
                          destination, 0.05)

    def isochrone(self, age=1 * 10 ** 8.0, metallicity=0.0):
        """
        This function adds several necessary
        columns (e.g. isWR, logg, and
        whether a star is WR) to the existing
        isochrone files' tables.
        Important note on convention: I will specify
        the phase of the star as follows:
        White dwarf -> 101
        All other stars (no neutron stars or black holes) -> -1
        Just to make sure, I do have word that neutron stars
        and black holes are NOT INCLUDED in the BPASS
        evolution grids.
        If you are REALLY curious about black holes, try using TUI
        with BPASS.

        Parameters
        ----------
        dir_in: string
        The string representation of the absolute path to the directory
        of preprocessed isochrones.This is for customization/Testing purposes.
        For folks do not want to store in the models directory.
        age: float or double
        The age of the isochrone in years
        metallicity: float or double
        The log10(Z/Z_solar) of the star
        (HERE THIS IS IN terms of Z which is the initial
        metallicity mass fraction from the formula introduced
        in Grevesse & Noels 1993.)
        Z_solar = 0.020
        Output
        ------
        A more informative version of the preprocessed isochrone.
        New Columns:
        isWR: Boolean
            whether a single or primary is a WR star
        isWR2: Boolean
        phase: Integer (101 if white dwarf, 12 otherwise)
              indicates whether a star is a white dwarf
        logg: float or double
             the log10(surface gravity of the primary/single in m/s^2)
        logg2: float or double
             the log10(surface gravity of the secondary in m/s^2)
             set to N.A. when the system is nonbinary
        log_ai: float or double
             the log10(separation of the primary and secondary in AU)
             set to N.A. if the system is not a binary system
        More clearly named are the columns for the
        current masses of both stars.
        'mass_current' column stands for the current
        mass of the primary star
        'mass_current2' column stands for the current
        mass of the secondary star
        """
        oldmetallicity = metallicity
        # The following metallicity fraction as how BPASS is organized.
        # Borrow from matt hosek's code:
        # Look for the closest valid age to the input age
        metallicity = self.z_solar * 10 ** metallicity
        log_age = math.log10(age)
        if log_age < np.min(self.age_list) or log_age > np.max(self.age_list):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
        if (metallicity < np.min(self.z_list) or
                metallicity > np.max(self.z_list)):
            logger.error('Requested metallicity {0} is out of bounds.'.
                         format(metallicity))
        iso_file = 'iso' + str(round(log_age, 1)) + '.fits'
        # Find the closest valid metallicity to the given metallicity
        closest_metallicity = min([x for x in self.z_file_map],
                                  key=lambda x: abs(metallicity - x))
        z_dir = self.z_file_map[closest_metallicity]
        # Use the full path to the desired isochrone file.

        full_iso_file = self.models_dir + 'iso/' + z_dir + iso_file
        iso = Table.read(full_iso_file, format='fits')
        # Create columns for isWR, logg,
        # Phase. I will need to find how
        # phase for individual stars will be added
        # then I will fill it in with appropriate values.
        isWR = Column(np.repeat(False, len(iso)), name='isWR')
        isWR2 = Column(np.repeat(False, len(iso)), name='isWR2')
        colG = Column(np.repeat(np.nan, len(iso)), name='logg')
        colP = Column(np.repeat(np.nan, len(iso)), name='phase')
        colP2 = Column(np.repeat(np.nan, len(iso)), name='phase2')
        collog_a = Column(np.repeat(np.nan, len(iso)), name='log_ai')
        iso.add_column(colG)
        iso.add_column(isWR)
        iso.add_column(colP)
        iso.add_column(colP2)
        iso.add_column(isWR2)
        iso.add_column(collog_a)
        # We may as well delete a for loop here
        # Convert to CGS in the next two lines.
        iso['logg'] = np.log10((iso['M1'] * cs.GM_sun /
                                (((10 ** iso['log(R1)']) *
                                  cs.R_sun) ** 2)) *
                               un.s * un.s/un.cm)
        iso['logg2'] = np.log10((iso['M2'] * cs.GM_sun /
                                (((10 ** iso['log(R2)']) *
                                  cs.R_sun) ** 2)) *
                                un.s * un.s/un.cm)
        # Apply Kepler's Third law to find the log of initial separation of the
        # binary system.
        iso['log_ai'] = np.log10((cs.GM_sun * (un.s ** 2) / (un.m ** 3)) **
                                 (1 / 3) *
                                 ((iso['mass'] + iso['mass2']) *
                                  ((24 * (60 ** 2) * (10 **
                                                      iso['initl_logP'])) **
                                   2) / (4 * (np.pi) ** 2)) ** (1 / 3) *
                                 (1 / cs.au) * un.m)
        iso.rename_column('M1', 'mass_current')
        iso['age'] = np.log10(iso['age'])

        # Using Stanway and Elridge Criterion for calculating whether
        # a star is a WR star or not for using the PotsDam Atmospheres
        # Decided to get rid of loops in order to take fuller advantage
        # of the C-based numpy!
        iso['isWR'] = (iso['X'] < 0.40) & (iso['log(T1)'] >= 4.45)
        iso['isWR2'] = (iso['X'] < 0.40) & (iso['log(T2)'] >= 4.45)
        iso['phase'] = 5
        iso['phase2'] = 5
        iso['phase'][np.where((iso['logg'] > 6.9) & (iso['log(L1)'] < -1) &
                              (iso['mass_current'] < 1.4))[0]] = 101
        iso['phase'][np.where(((iso['source'] == 2) | (iso['source']==3) |
                               (iso['source'] == 4)) &
                              (np.round(iso['mass_current'], 1) == 1.4))[0]] = 102
        iso['phase'][np.where(((iso['source'] == 2) | (iso['source'] == 3) |
                               (iso['source'] == 4)) &
                              (np.round(iso['mass_current']) > 3.0))[0]] = 103
        iso['phase2'][np.where((iso['logg2'] > 6.9) &
                               (iso['log(L2)'] < -1) &
                               (iso['M2'] < 1.4))[0]] = 101
        # Changing column name to
        # Making sure that names of the columns are consistent with
        # general format of isochrone.
        iso.rename_column('log(T1)', 'logT')
        iso.rename_column('M2', 'mass_current2')
        iso.rename_column('log(L1)', 'logL')
        iso.rename_column('age', 'logAge')
        iso['logT'][np.where(iso['phase'] == 103)] = -np.inf
        iso['logL'][np.where(iso['phase'] == 103)] = -np.inf
        iso.meta['log_age'] = age
        iso.meta['metallicity_in'] = oldmetallicity
        iso.meta['metallicity_act'] = np.log10(closest_metallicity /
                                               self.z_solar)
        return iso

# === Code from IsochroneMakerReformattedVersionNew === #
# Made to make Pytest shut the heck up about modules not
# being able to be found

import pandas as pd
from astropy.io import fits
from astropy.table import Table
import os
import hoki
from hoki import load
import glob

possible_secondary_q = ['0.1', '0.2', '0.3',
                        '0.4', '0.5',
                        '0.6', '0.7', '0.8', '0.9']
# The list encompassing possible primary
# mass for primary and secondary. Trying to encompass
# all possible cases, even if there is a bit of redundancy and excess.
mass_list = [round(0.1, 1)] + [round(0.12 + x * 0.020, 2) for x in range(95)]
mass_list = mass_list + [round(2.05 + x * 0.05, 2) for x in range(20)]
mass_list = mass_list + [round(3.1 + 0.1 * x, 1) for x in range(70)]
mass_list = (mass_list + [11 + x for x in range(90)] +
             [120, 400, 500] + [125 + 25 * x for x in range(8)])


def convert(x):
    """
    Helper function to use for possible mass list. If the
    number can be represented by an integer,
    it turns into an integer.
    This is to make sure that I am really matching masses.
    """
    if (x % 1 == 0.0):
        return int(x)
    else:
        return x


mass_list = list(map(convert, mass_list))
period_list = ['0', '0.2', '0.4', '0.6', '0.8', '1', '1.2', '1.4',
               '1.6', '1.8', '2', '2.2', '2.4', '2.6', '2.8', '3',
               '3.2', '3.4', '3.6', '3.8', '4']
combos = []
for y in period_list:
    for z in possible_secondary_q:
        combos.append((z, y))


def assign_props(dictionary, input_str, y):
    dictionary[input_str] = y
    return input_str


# Source: https://stackoverflow.com/questions/44369504/
# how-to-convert-entire-dataframe-values-to-float-in-pandas
vals = hoki.dummy_dict.values()
vals = [hoki.dummy_dict[keyword] for keyword in ['timestep', 'age',
                                                 'log(R1)', 'log(T1)',
                                                 'log(L1)', 'M1', 'X',
                                                 'P_bin', 'log(a)', 'M2',
                                                 'log(R2)', 'log(T2)',
                                                 'log(L2)']]
cols_to_keep = ["col" + str(v + 1) for v in vals]
# According to the BPASS v2.2.1 manual, much of
# columns 50 and onward
# are basically spectra and atmosphere related models.
# They will be calculated later.
# I have only kept the column numbers corresponding to
# Source for inverse mapping: https://stackoverflow.com/questions/483666/
# reverse-invert-a-dictionary-mapping
# I create a mapping from column NUMBER to
# column name
invmap = {u: v for v, u in hoki.dummy_dict.items()}
# Accounting for zero indexing, write out a list
# of column names to assign to each BPASS dat file
# column name.
# Keep in mind the non-zero indexing of the data
# table itself.
lisnp = [invmap[int(x[3:]) - 1] for x in cols_to_keep]


def find_mergers(star_table, initial_mass):
    """
    given star_table, an Astropy Table of primary masses of a
    NEWBINMODS stellar model
    and age, and given the initial mass of the stellar model,
    sorts the star_table by age and
    returns the index at which the primary star's mass increases.
    If the mass of the primary star ends up being greater than the
    initial mass at the very first
    age represented by the star_table, then 0 is designated as
     the index where the merger occurs.
    If the merger does not occur then the length of the star_table
     is returned. (cannot be a valid index
    in a table). This is to ensure that all rows corresponding
     to positions after the returned index of the
    star_table are merger.
    """

    if (star_table['M1'][0] > initial_mass + 0.01):
        return 0
    # returns index where the merger is completed
    # assuming that f is a merger model.
    # We need to make sure that Mass of primary is increasing as age increases.
    # I make sure that age is actually increasing and
    # that the mass of the primary
    # is increasing by at least 0.01 solar masses.
    # Per Dr. Jan Edlridge, the increase in mass is supposed to be
    # an indicator of a merger happening
    val = np.where((np.diff(star_table['M1']) > 0.01) &
                   (np.diff(star_table['age']) >= 0))[0]
    if (len(val)):
        return val[-1] + 1
    return len(star_table)


def reformatter(destination, metallicity):
    """I will expedite the process of writing and saving files by
    making the BPASS stellar input files binary fits tables.
    This function reads every stellar model of the specified metallicities
    as an astropy table. Deletes photometry related columns of each stellar
    model and renames the columns to a human readable name as prescribed by
    the hoki.dummy_dict. Saves the reformatted BPASS stellar evolution model
    as a FITS formatted file for the sake of speed of processing later on.

    Parameters
    ----------
    destination: String
        describes the absolute path to the directory of BPASS models.
        Has a dash at the end
    metallicity: List of strings
        represents the metallicities of the BPASS input files that are
        to be processed.

    Remarks: Run this on Terminal not on Ipython. This function goes
    through tens of thousands of files full of data; expect it to take
    a long time to complete during your setup if you choose to use this
    via Jupyter Notebook. Using the reformatter on the the
    terminal will help save time.
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
    # I will use the Python HashSet since (assuming that strings
    # have a good HashCode, adding to the HashSet should be fast in
    # normal case.
    if (not os.path.isdir(destination)):
        os.makedirs(destination)
    setofAll = set()
    for x in metallicity:
        # Generate all possible BPASS models with the given metallicities.
        # Get a set of all single star models.
        a = set(glob.glob("{}/NEWSINMODS/{}/*".format(hoki.MODELS_PATH, x)))
        # Get a set of all Quasi-Chemicall Homogeneous Evolution (Binary)
        # Models.
        d = set(glob.glob("{}/NEWSINMODS/{}hmg/*".format(hoki.MODELS_PATH, x)))
        # Get a set of all regular binary and merged binary models.
        b = set(glob.glob(("{}/NEWBINMODS/" +
                           "NEWBINMODS/{}/*").format(hoki.MODELS_PATH, x)))
        # Get a set of all models where the primary is a compact remnant.
        c = set(glob.glob(("{}/NEWBINMODS/" +
                          "NEWSECMODS/{}_2/*").format(hoki.MODELS_PATH, x)))
        setofAll = setofAll.union(a, b, c, d)
        if not os.path.isdir("{}/{}/".format(destination, x)):
            os.makedirs("{}/{}/".format(destination, x))
        # Creates container directory for stellar models grouped by metallicity
        # and whether they are single or binary related
        # the latter included mergers, QHE Models
        # Single star models that are not the hmg models
        # go into the <Metallicity>/sin
        # directory.
    new_sin_to_met = len(hoki.MODELS_PATH) + len("/NEWSINMODS/")
    new_sin_and_met = len(hoki.MODELS_PATH) + len("/NEWSINMODS/zxxx/")
    new_hmg_and_met = len(hoki.MODELS_PATH) + len("/NEWSINMODS/zxxxhmg/")
    new_bin_to_met = len(hoki.MODELS_PATH) + len("/NEBINMODS/NEWBINMODS/")
    new_sec_to_met = len(hoki.MODELS_PATH) + len("/NEBINMODS/NEWSECMODS/")
    new_bin_and_met = len(hoki.MODELS_PATH) + \
    len("/NEBINMODS/NEWBINMODS/zxxx/")
    new_sec_and_met = len(hoki.MODELS_PATH) + \
    len("/NEBINMODS/NEWBINMODS/zxxx_2/")
    for x in setofAll:
        astroTable = Table.read(x, format='ascii')
        # Drop all of the columns that are not defined by the
        # HOKI package provided mapping of column name to
        # column number.
        # We will be recalculating the photometry related stuff
        # Using zero indexing for list of columns,
        # All columns (except column number 12, 13, 36)
        # up to 48 inclusive will be kept
        # when we create the isochrone.
        astroTable.keep_columns(cols_to_keep)
        astroTable.rename_columns(astroTable.colnames, lisnp)
        astroTable.sort("age")
        # Will remove photometry related columns
        # Merger Related indicates whether the star is
        # a model with a merger in it.
        # Checking if the model is a single star model.
        sin_segment = x[new_sin_to_met:new_sin_to_met + 4]
        # SEE IF THE FILE IS IN THE NEWSINMODS directory
        if (x[len(hoki.MODELS_PATH) + 1:len(hoki.MODELS_PATH) +
                len("NEWSINMODS") + 1] == "NEWSINMODS" and x[new_sin_to_met +
                                                         4:new_sin_to_met +
                                                         7] != "hmg"):

            astroTable.write(destination +
                             "/" + sin_segment + "/Fits" +
                             "Models{}sin.fits".format(
                                    x[new_sin_and_met:]),
                             format='fits', overwrite=True)
        else:
            # Checking if the model is a HMG Binary model
            # (Still for some reason in SINMODS)
            if x[len(hoki.MODELS_PATH) + 1:len(hoki.MODELS_PATH) + 1 +
                 len("NEWSINMODS")] == "NEWSINMODS":

                astroTable.write(destination + "/" +
                                 sin_segment + "/Fits" +
                                 "Models{}hmg.fits".format(x[new_hmg_and_met:]),
                                 format='fits', overwrite=True)
                continue
                # Models that count as secondary star with compact remnants
                # go into the <metallicity>/bin subdirectory of
                # the destination.
            if x[len(hoki.MODELS_PATH) + 1:len(hoki.MODELS_PATH) + 1 +
                 len("NEWBINMODS/XXXXXXXXXX")] == "NEWBINMODS/NEWSECMODS":

                astroTable.write(destination +
                                 "/" + x[new_sec_to_met:
                                         (new_sec_and_met) - 2] +
                                 "/FitsModels{}sec.fits".format(x[new_sec_and_met + 1:]),
                                 format='fits', overwrite=True)
            else:
                # Models that count as secondary star with compact remnants
                # go into the <metallicity>/ subdirectory of the destinatuion.

                astroTable.write(destination + '/' +
                                 x[new_bin_to_met: new_bin_and_met] +
                                 "/FitsModels{}bin.fits".format(x[new_bin_and_met + 1:]),
                                 format='fits', overwrite=True)


def extractor(age, metallicity, input_dir, bpass_evo_dir,
              margin):
    """
    The following extracts the preprocessed BPASS
    isochrone file
    from the pool of listed files. From the <metallicity> subdirectory
    of the input directory, the function inspects all of the files
    and sees if there are rows in the files that correspond to
    the specified age (log age to be precise).
    Checks for merger related (models_type=0) models
    and mark them using the functions merger_related and
    find_mergers respectively.



    parameters
    ----------
    age: float or double
          The age in years of the isochrone desired
    metallicity: string
          The BPASS metallicity of the isochrone desired
    input_dir: string
          The absolute path to the directory with the reformatted BPASS stellar
          models. The path should be to a folder whose name was recently
          passed into reformatter.
    bpass_evo_dir: string
          The absolute path to the directory where the outputs of this program,
          preprocessed isochrone files for BPASS, will
          be stored in subdirectories
          named after their metallicities.
    margin: float or double
             Contains the margin of error of base-10 logAge that will be
             tolerated for rows of stellar model that are included in
             the isochrone.

    Output:
    Saves into the specified destination
    directory a reformatted, FITS format version
    of a preprocessed BPASS isochrone of the specified age
    and metallicity whose individual
    rows correspond to stars of age closest to
    the specified age and within the specified
    margin of error of log(age of the star in years)

    Parameters added to the table:
    ‘mass’ - The column corresponding to the initial mass
             of the primary/single star.
    ‘mass2’ -  The column corresponding to the initial mass of
               the secondary. Set to nan if the system is not a binary
    ‘log_P’ - The column corresponding to the log(Period in days) of
              the binary. Set to nan if the system is not a binary
    ‘Single’ - The column that tells whether the stellar system is single
    ‘Merged?’ - Whether the binary star has been merged.

    """
    # Set of files that exist and have been caught by the looping.
    caught_no = set()
    # Mapping of file name to a tuple of properties
    # (primary star mass,(secondary star mass, initial logP for bin.fits types)
    names_to_prop = {}
    # Find all filenames of reformatted BPASS models from reformatter
    for x in mass_list:
        # Find all NEWBINMODS systems of the specified metallicity
        for y in combos:
            st = "{}/{}/FitsModelssneplot-{}-{}-{}-{}bin.fits"
            st = st.format(input_dir,metallicity,
                           metallicity,
                           str(x), y[0],
                           y[1])
            names_to_prop[st] = (str(x), y[0], y[1])
            if (os.path.isfile(st)):
                caught_no.add(st)
        # Find all single star systems of the specified metallicity
        st = "{}/{}/FitsModelssneplot-{}-{}sin.fits"
        st = st.format(input_dir, metallicity,
                       metallicity, str(x))
        names_to_prop[st] = (str(x), )
        if (os.path.isfile(st)):
            caught_no.add(st)
        # Find all Binary QHE systems of the specified metallicity
        st = "{}/{}/FitsModelssneplot-{}-{}hmg.fits"
        st = st.format(input_dir, metallicity,
                       metallicity, str(x))
        names_to_prop[st] = (str(x), )
        if (os.path.isfile(st)):
            caught_no.add(st)
        # Find all NEWSINMODS systems of the specified metallicity
        st = "{}/{}/FitsModelssneplot_2-{}-{}-*sec.fits"
        st = st.format(input_dir, metallicity,
                       metallicity,
                       str(x))
        li = glob.glob(str(st))
        sec_files = set(li)
        caught_no = caught_no.union(sec_files)
        [assign_props(names_to_prop, name, (x,)) for name in sec_files]
    len_of_heading = len("{}/{}/FitsModelssneplot_2-{}-".format(input_dir, metallicity, metallicity))
    bigOne = None
    initial = True
    initlMass = np.nan
    initlMass2 = np.nan
    indicesOfInterest = None
    entries = glob.glob(str("{}/{}/*".format(input_dir, metallicity)))
    suffix_len = len("xxx.fits")  # 8
    for x in entries:
        # If I want to see the name of the file, uncomment below.
        # print(x)
        # x is the name of the reformatted stellar evolution file.
        # Rest carries the initial mass of the system and for NEWBINMODS
        # systems carries the secondary and log_P in days
        rest = names_to_prop[x]
        # The following conditional would have been that x is a path to a file,
        # but this seems to be always true.
        # as I am getting the entries from globs and
        # prior checks as to whether they are existing files.
        # We only want to consider the
        org = Table.read(x, format='fits')
        indicesOfInterest = np.where(np.abs(np.log10(org['age']) -
                                            age) <= margin)[0]
        f = org[indicesOfInterest]  # f stands for frame in DataFrame
        indicesOfInterest = np.array(indicesOfInterest)
        # If no stars have a log10(age) within margin of
        # the given lage in log-10 years, we
        # must skip over to the next model.
        if (len(f) != 0):
            # Find the star with age closest to the input log(Age of the star)
            filterDown = np.where(np.abs(f['age'] - 10 ** age) ==
                                  np.min(np.abs(f['age'] - 10 ** age)))[0]
            f = f[filterDown]
            indicesOfInterest = np.array([indicesOfInterest[filterDown][0]])
            if len(f) != 0:
                f = f[[0]]
                f['source'] = 0
                indexlen = 1
                if (x[-suffix_len:-5] == 'bin'):
                    f['single'] = np.repeat(False, indexlen)
                    initlMass = float(rest[0])
                    f['mass'] = np.repeat(initlMass, indexlen)
                    f['mass2'] = np.repeat(initlMass * float(rest[1]),
                                           indexlen)
                    f['initl_logP'] = np.repeat(float(rest[2]), indexlen)
                    # Now, for binaries, I check whether a model is
                    # a merger model
                    # is the same for all rows of the model
                    merge_pt = find_mergers(org[['age', 'M1']], initlMass)
                    # Find whether we should treat our model as
                    # one big star or still two stars
                    # will still keep initial parameters
                    # for the sake of working
                    # with the Duchene-Krauss distributions
                    f['mergered?'] = indicesOfInterest >= merge_pt
                    f['source'] = 1
                elif (x[-suffix_len:-5] == 'sec'):
                    # Recall that the ending of the file
                    # is going to be XXX.fits
                    # The XXX can be sec, hmg, bin, sin.
                    # I will explot the pattern that the log_P is
                    # going to be either 9 characters,
                    # 8 characters, or 7 characters
                    # long.
                    f['single'] = np.repeat(False, indexlen)
                    initlMass = float(rest[0])
                    # See if the number would begin
                    # right after a dash in index -18.
                    # 8 characters given for the xxx.fits
                    # 9 spaces for the decimal number.
                    # We want the dash before that as our cue.
                    if (x[-10 - 1 * suffix_len] == "-"):
                        # Here the decimal is number is 9
                        # characters long and
                        # xxx.fits has length of 8
                        log_P_in_days = float(x[-9 - suffix_len:
                                                -1 * suffix_len])
                        # Accounting for the dashes, find the
                        # sandwiched initial
                        # mass of the "remnant primary"
                        # sandwiched between initial secondary star mass
                        # and the initial logP
                        initlMass2 = float(x[(len_of_heading +
                                              len(str(rest[0]) + "-")):
                                             -10 - 1 * suffix_len])
                        f['source'] = 2
                    else:
                        # See if the number begins right
                        # after a dash in index -17
                        # 8 for the xxx.fits and 8 for
                        # the log_P and we want
                        # the preceding dash.
                              
                        if (x[-9 - 1 * suffix_len] == "-"):
                            
                            # Here the decimal is number is 8
                            # characters long and
                            # xxx.fits has length of 8.
                            log_P_in_days = float(x[-8 - suffix_len:
                                                    -1 * suffix_len])
                            # Accounting for the dashes, find
                            # the sandwiched initial mass of the
                            # "remnant primary"
                            # sandwiched between initial secondary star
                            # mass and the initial logP
                            initlMass2 = float(x[(len_of_heading +
                                                  len(str(rest[0]) +
                                                      "-")):
                                                 -9 - 1 * suffix_len])
                            # Assume that the number would begin
                            # at index -16. (From my observation, we only
                            # have a few choices for number of digits of
                            # the log_period
                            f['source'] = 3
                        else:
                            # Here the decimal is number is 7
                            # characters long and
                            # xxx.fits has length of 8
                            log_P_in_days = float(x[-7 - suffix_len: 
                                                    -1 * suffix_len])
                            # Accounting for the dashes on both
                            # ends of the "remnant primary" mass, find the
                            # sandwiched initial mass of the remnant primary
                            # sandwiched between initial secondary star mass
                            # and the initial logP
                            initlMass2 = float(x[(len_of_heading +
                                                  len(str(rest[0]) +
                                                      "-")):
                                                 -8 - 1 * suffix_len])
                            f['source'] = 4
                    f['mass'] = np.repeat(initlMass2, indexlen)
                    f['mass2'] = np.repeat(initlMass, indexlen)
                    temp1 = f['M2']
                    temp2 = f['log(T2)']
                    temp3 = f['log(R2)']
                    temp4 = f['log(L2)']
                    f['M2'] = f['M1']
                    f['log(R2)'] = f['log(R1)']
                    f['log(T2)'] = f['log(T1)']
                    f['log(L2)'] = f['log(L1)']
                    f['M1'] = temp1
                    f['log(R1)'] = temp3
                    f['log(T1)'] = temp2
                    f['log(L1)'] = temp4
                    f['initl_logP'] = np.repeat(log_P_in_days, indexlen)
                    f['mergered?'] = np.repeat(False, indexlen)
                    f['single'] = np.repeat(False, indexlen)
                elif (x[-1 * suffix_len: -5] == 'hmg'):
                    f['single'] = np.repeat(True, indexlen)
                    initlMass = float(rest[0])
                    # To be consistent with obtaining initial mass
                    # (or whichever value is closest) for the HMG Model).
                    # Per Eldridge the QHE models are single star models.
                    initlMass2 = np.nan
                    P_in_days = np.nan
                    f['mass'] = np.repeat(initlMass, indexlen)
                    f['mass2'] = np.repeat(initlMass2, indexlen)
                    f['initl_logP'] = np.repeat(P_in_days, indexlen)
                    f['mergered?'] = np.repeat(False, indexlen)
                    f['source'] = 5
                else:
                    f['single'] = np.repeat(True, indexlen)
                    # Single stars do not have companions.
                    initlMass = rest[0]
                    initlMass2 = np.nan
                    P_in_days = np.nan
                    f['mass'] = np.repeat(initlMass, indexlen)
                    f['mass2'] = np.repeat(initlMass2, indexlen)
                    f['initl_logP'] = np.repeat(np.nan, indexlen)
                    f['mergered?'] = np.repeat(False, indexlen)
                    f['source'] = 6
                if initial:
                    initial = False
                    bigOne = f.to_pandas()
                else:
                    bigOne = pd.concat([f.to_pandas(), bigOne])
    if not isinstance(bigOne, type(None)) and not (bigOne['age'].empty):
        bigOne = bigOne.apply(pd.to_numeric, errors='coerce')
        reduced = Table.from_pandas(bigOne)
        if not os.path.isdir('{}/iso/{}/'.format(bpass_evo_dir, metallicity)):
            os.makedirs('{}/iso/{}/'.format(bpass_evo_dir, metallicity))
        reduced.write("{}iso/{}/iso{}.fits".format(bpass_evo_dir,
                                                   metallicity, str(age)),
                      format='fits',
                      overwrite=True)
