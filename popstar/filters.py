import numpy as np
from astropy.table import Table
import pysynphot
import warnings
import os
import pdb

# Set path to filter functions
code_dir = os.path.dirname(__file__)
filters_dir = code_dir[:-8]+'/filt_func/'

def get_nirc2_filt(name):
    """
    Define nirc2 filter as a pysynphot spectrum object
    """
    # Read in filter info
    try:
        t = Table.read('{0}/nirc2/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find NIRC2 filter file {0}/nirc2/{1}.dat'.format(filters_dir, name))

    wavelength = t[t.keys()[0]]
    transmission = t[t.keys()[1]]

    # Lets fix wavelength array for duplicate values
    diff = np.diff(wavelength)
    idx = np.where(diff <= 0)[0]

    while len(idx) != 0:
        wavelength[idx+1] += 1.0e-8
        
        diff = np.diff(wavelength)
        idx = np.where(diff <= 0)[0]
        print( 'Duplicate entry loop' )

    # Get rid of all entries with negative transmission
    idx = np.where(transmission > 1)[0]

    # Convert wavelength to Angstroms, transmission to ratio
    wavelength = wavelength[idx] * 10**4
    transmission = transmission[idx] / 100.0 # convert from % to ratio

    # Make spectrum object
    spectrum = pysynphot.ArrayBandpass(wavelength, transmission, waveunits='angstrom',
                                           name='NIRC2_{0}'.format(name))

    return spectrum

def get_2mass_filt(name):
    """
    Define the 2mass filters as a pysynphot spectrum object
    """
    # Read in filter info
    try:
        t = Table.read('{0}/2mass/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find 2MASS filter file {0}/2mass/{1}.dat'.format(filters_dir, name))

    wavelength = t[t.keys()[0]]
    transmission = t[t.keys()[1]]

    # Convert wavelength to Angstroms
    wavelength = wavelength * 10**4

    # Make spectrum object
    spectrum = pysynphot.ArrayBandpass(wavelength, transmission, waveunits='angstrom',
                                           name='2MASS_{0}'.format(name))

    return spectrum
    

def get_vista_filt(name):
    """
    Define vista filter as pysynphot spectrum object
    """
    # Read in filter info
    try:
        t = Table.read('{0}/vista/VISTA_Filters_at80K_forETC_{1}.dat'.format(filters_dir, name),
                           format='ascii')
    except:
        raise ValueError('Could not find VISTA filter file {0}/vista/VISTA_Filters_at80K_forETC_{1}.dat'.format(filters_dir, name))    

   # Wavelength must be in angstroms, transmission in fraction
    wave = t['col1'] * 10
    trans = t['col2'] * 0.01
    
    # Change any negative numbers to 0, as well as anything shortward
    # of 0.4 microns or longward of 2.9 microns
    # (no VISTA filter transmissions beyond these boundaries)
    bad = np.where( (trans < 0) | (wave < 4000) | (wave > 29000) )
    trans[bad] = 0
    
    # Now we can define the VISTA filter bandpass objects
    spectrum = pysynphot.ArrayBandpass(wave, trans, waveunits='angstrom', name='VISTA_{0}'.format(name))
    
    return spectrum

def get_decam_filt(name):
    """
    Define DECAM filter as pysynphot object
    """
    # Read in filter info
    try:
        t = Table.read('{0}/decam/DECam_filters.txt'.format(filters_dir), format='ascii')
        t.rename_column('Y', 'y')
        
        cols = np.array(t.keys())
        idx = np.where(cols == name)[0][0]

        trans = t[cols[idx]]
    except:
        raise ValueError('Could not find DECAM filter {0} in {1}/decam/DECam_filters.txt'.format(name, filters_dir))         

    # Limit to unmasked regions only
    mask = np.ma.getmask(trans)
    good = np.where(mask == False)
    
    # Convert wavelengths from nm to angstroms, while eliminating masked regions
    wave = t['wavelength'][good] * 10.
    trans = trans[good]
    wave = np.ma.filled(wave)
    trans = np.ma.filled(trans)

    spectrum = pysynphot.ArrayBandpass(wave, trans, waveunits='angstrom', name='decam_{0}'.format(name))

    return spectrum

def get_PS1_filt(name):  
    """
    Define PS1 filter as pysynphot object
    """
    try:
        t = Table.read('{0}/ps1/PS1_filters.txt'.format(filters_dir), format='ascii')
        t.rename_column('col1', 'wave')
        t.rename_column('col2', 'open')
        t.rename_column('col3', 'g')
        t.rename_column('col4', 'r')
        t.rename_column('col5', 'i')
        t.rename_column('col6', 'z')
        t.rename_column('col7', 'y')

        cols = np.array(t.keys())
        idx = np.where(cols == name)[0][0]

        trans = t[cols[idx]]
    except:
        raise ValueError('Could not find PS1 filter {0} in {1}/ps1'.format(name, filters_dir))         

    # Convert wavelengths from nm to angstroms
    wave = t['wave'] * 10.
    
    spectrum = pysynphot.ArrayBandpass(wave, trans, waveunits='angstrom', name='ps1_{0}'.format(name))

    return spectrum

def get_jwst_filt(name):
    """
    Define JWST filter as pysynphot object
    """
    try:
        t = Table.read('{0}/jwst/{1}.txt'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find JWST filter {0} in {1}/jwst'.format(name, filters_dir))         

    # Convert wavelengths to angstroms
    wave = t['microns'] * 10**4.
    trans = t['throughput']

    # Change any negative numbers to 0
    bad = np.where(trans < 0)
    trans[bad] = 0
    
    spectrum = pysynphot.ArrayBandpass(wave, trans, waveunits='angstrom', name='jwst_{0}'.format(name))

    return spectrum    

def get_Johnson_Glass_filt(name):
    """
    Define Johnson-Glass filters as pysynphot object
    """
    try:
        t = Table.read('{0}/Johnson_Glass/{1}.txt'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find Johnson-Glass filter {0} in {1}/Johnson_Glass'.format(name, filters_dir))         

    # Convert wavelengths to angstroms
    wave = t['col1'] * 10.
    trans = t['col2']

    # Change any negative numbers to 0
    bad = np.where(trans < 0)
    trans[bad] = 0
    
    spectrum = pysynphot.ArrayBandpass(wave, trans, waveunits='angstrom', name='jg_{0}'.format(name))

    return spectrum    

def get_nirc1_filt(name):
    """
    Define Keck/NIRC filters as pysynphot object
    """
    try:
        t = Table.read('{0}/nirc1/{1}.txt'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find NIRC1 filter {0} in {1}/nirc1'.format(name, filters_dir))         

    # Convert wavelengths to angstroms
    wave = t['col1'] * 10**4
    trans = t['col2']
    
    # Lets fix wavelength array for duplicate values or negative vals;
    # delete these entries
    diff = np.diff(wave)
    idx = np.where(diff <= 0)[0]

    while(len(idx) != 0):
        bad = idx + 1

        wave = np.delete(wave, bad)
        trans = np.delete(trans, bad)

        diff = np.diff(wave)
        idx = np.where(diff <= 0)[0]
        
    # Change any negative transmission vals to 0
    bad = np.where(trans < 0)
    trans[bad] = 0

    spectrum = pysynphot.ArrayBandpass(wave, trans, waveunits='angstrom', name='nirc1_{0}'.format(name))

    return spectrum    

def get_ctio_osiris_filt(name):
    """
    Define CTIO/OSIRIS filters as pysynphot object
    """
    try:
        t = Table.read('{0}/CTIO_OSIRIS/{1}.txt'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find CTIO/OSIRIS filter {0} in {1}/CTIO_OSIRIS'.format(name, filters_dir))         

    # Convert wavelengths to angstroms
    wave = t['col1'] * 10**4
    trans = t['col2']

    # Change any negative numbers to 0
    bad = np.where(trans < 0)
    trans[bad] = 0
    
    spectrum = pysynphot.ArrayBandpass(wave, trans, waveunits='angstrom', name='ctio_osiris_{0}'.format(name))

    return spectrum    
