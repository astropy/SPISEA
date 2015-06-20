#
#
#
import math
import logging
from numpy import searchsorted, genfromtxt
import numpy as np
import os
import glob

log = logging.getLogger('starEvo')

class StellarEvolution():
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
    
class GenevaStellarEvolution(StellarEvolution):
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

        super().__init__(model_dir, age_list, mass_list, z_list)

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
        iso_file = 'iso_' + str(self.age_list[age_idx]) + '.dat'
        
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

class EkstromStellarEvolution(StellarEvolution):
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
        #iso_file = 'iso_' + str(self.age_list[age_idx]) + '.dat'
        
        # find closest metallicity value
        #z_idx = searchsorted(self.z_list, z_defined, side='left')
        #z_dir = self.z_file_map[self.z_list[z_idx]]
        
        # generate isochrone file string
        #full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file
        
        # return isochrone data
        return genfromtxt(full_iso_file, comments='#')

    def format_isochrones(input_iso_dir):
        r"""
        Parse iso.dat (filename hardcoded) file downloaded from Ekstrom+12
        models, create individual isochrone files for the different ages.

        input_iso_directory should lead to Ekstrom2012/iso/<metallicity> directory,
        where iso.dat file should be located

        Creates two new directories, rot and norot, which contain their respective
        isochrones.
        """
        # Move into metallicity direcotry, read iso.dat file
        os.chdir(input_iso_dir)
    
        print 'Read Input: this is slow'
        iso = Table.read('iso.dat', format='ascii')
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
            tmp_r.write('rot/iso_{0:4.2f}.dat'.format(age), format='ascii')
            tmp_n.write('norot/iso_{0:4.2f}.dat'.format(age), format='ascii')

        # Return to some directory
        #os.chdir(<somewhere>)
        return

#---------------------------------------#
# Now for the Parsec version 1.2s models
#---------------------------------------#

class ParsecStellarEvolution(StellarEvolution):
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
        #iso_file = 'iso_' + str(self.age_list[age_idx]) + '.dat'
        
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
        os.chdir(input_iso_dir)
        
        # Work on each metallicity isochrones individually
        for metal in metallicity_list:
            # More into metallicity directory, read isochrone file
            os.chdir(metal)

            isoFile = glob.glob('output*')
            print 'Read Input: this is slow'
            iso = Table.read(isoFile[0], format='ascii')
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
                tmp.write('iso_{0:4.2f}.dat'.format(age), format='ascii')

            # Move back into iso directory
            os.chdir('..')

        # Return to starting directory
        #os.chdir(<somewhere>)
        return


#---------------------------------------#
# Now for the Pisa (Tognelli+11) models
#---------------------------------------#

class PisaStellarEvolution(StellarEvolution):
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
        #iso_file = 'iso_' + str(self.age_list[age_idx]) + '.dat'
        
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
        os.chdir(input_iso_dir)
        # Work on each metallicity directory individually
        for metal in metallicity_list:
            # Move into directory, check to see if files are already formatted
            os.chdir(metal)

            if os.path.exists('iso_6.00.dat'):
                print 'Files in {0:s} already formatted'.format(metal)
            else:
                # Create a ReadMe with the original file names to preserve the
                # model details
        
                cmd = "ls *.DAT > ReadMe"
                os.system(cmd)

                # Collect all filenames in a list, rename files one
                # by one
                isoFile_list = glob.glob('*.DAT')
                for File in isoFile_list:
                    name = File.split('_')
                    # Extract iso age from filename
                    age = float(name[1][1:])
                    logAge = np.log10(age * 10**6)

                    cmd = "mv {0:s} iso_{1:4.2f}.dat".format(File, logAge)
                    os.system(cmd)

            # Return to overhead directory
            os.chdir('..')
        # Return to some direction
        #os.chdir(<somewhere>)
        return
