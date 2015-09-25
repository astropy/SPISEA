import logging
import numpy as np
import pysynphot
import os
import glob
from astropy.io import fits
from astropy.table import Table
import time
import pdb

log = logging.getLogger('atmospheres')

def get_kurucz_atmosphere(metallicity=0, temperature=20000, gravity=4):
    """
    metallicity in [Fe/H] (def = +0.0): 
          +1.0, +0.5, +0.3, +0.2, +0.1, +0.0, -0.1, -0.2, -0.3, -0.5, -1.0, 
          -1.5, -2.0, -2.5, -3.0, -3.5, -4.0, -4.5, -5.0.

    temperatures (def = 20,000 Kelvin):
                Temperature Range      Grid Step
                       K                   K

                  3000 - 10000            250 
                 10000 - 13000            500
                 13000 - 35000           1000
                 35000 - 50000           2500

    log gravity (def = 0.0) in the range of 0.0 - 5.0 in 0.5 increments
    """
    sp = pysynphot.Icat('k93models', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp.flux != 0)[0]
    if len(idx) == 0:
        print 'Could not find Kurucz 1993 atmosphere model for'
        print '  temperature = %d' % temperature
        print '  metallicity = %.1f' % metallicity
        print '  log gravity = %.1f' % gravity

    return sp

def get_castelli_atmosphere(metallicity=0, temperature=20000, gravity=4):
    """
    metallicity in [Fe/H] (def = +0.0): 
          0.0, -0.5, -1.0, -1.5, -2.0, -2.5.

    temperatures (def = 20,000 Kelvin):
                Temperature Range      Grid Step
                       K                   K

                  3000 -  13000            250 
                  13000 - 50000            1000

    log gravity (def = 4.0) in the range of 2.0 - 5.0 in 0.5 increments
    """
    # Round given temp, gravity to closest model that exists
    # No models below 3000 K
    if temperature < 3000:
        print 'No ck04 model below 3000 K'

    # Catch the edge cases for the hotest stars where the newest 
    # evolution models have logg < 4.0 but the atmosphere models
    # aren't available. HACK! FIX!
    if (temperature > 50000):
        print 'Changing temp from {0:5.0f} to 50000'.format(temperature)
        temperature = 50000
    if (temperature > 49000) and (gravity < 5.0):
        print 'Changing gravity for T=', temperature, ' logg=', gravity
        gravity = 5.0
    if (temperature > 39000) and (gravity < 4.5):
        print 'Changing gravity for T=', temperature, ' logg=', gravity
        gravity = 4.5
    if (temperature > 31000) and (gravity < 4.0):
        print 'Changing gravity for T=', temperature, ' logg=', gravity
        gravity = 4.0
    if (temperature > 26000) and (gravity < 3.5):
        print 'Changing gravity for T=', temperature, ' logg=', gravity
        gravity = 3.5
    if (temperature > 19000) and (gravity < 3.0):
        print 'Changing gravity for T=', temperature, ' logg=', gravity
        gravity = 3.0
    if (temperature > 11750) and (gravity < 2.5):
        print 'Changing gravity for T=', temperature, ' logg=', gravity
        gravity = 2.5

    if (temperature > 9000) and (gravity < 2.0):
        print 'Changing gravity for T=', temperature, ' logg=', gravity
        gravity = 2.0
    if (temperature > 8250) and (gravity < 1.5):
        print 'Changing gravity for T=', temperature, ' logg=', gravity
        gravity = 1.5
    if (temperature > 7500) and (gravity < 1.0):
        print 'Changing gravity for T=', temperature, ' logg=', gravity
        gravity = 1.0

    # Also edge case where gravity > 5.0, set to gravity = 5.0. This
    # is true at all temperatures. HACK!
    if gravity > 5.0:
        print 'Changing gravity for T=', temperature, ' logg=', gravity
        gravity = 5.0
    if (gravity < 0.5):
        print 'Changing gravity for T=', temperature, ' logg=', gravity
        gravity = 0.5
        
        
    sp = pysynphot.Icat('ck04models', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp.flux != 0)[0]
    if len(idx) == 0:
        print 'Could not find Castelli and Kurucz 2004 atmosphere model for'
        print '  temperature = %d' % temperature
        print '  metallicity = %.1f' % metallicity
        print '  log gravity = %.1f' % gravity

    return sp

def get_nextgen_atmosphere(metallicity=0, temperature=5000, gravity=4):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 5000)
    gravity = log gravity (def = 4.0)
    """
    sp = pysynphot.Icat('nextgen', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp.flux != 0)[0]
    if len(idx) == 0:
        print 'Could not find NextGen atmosphere model for'
        print '  temperature = %d' % temperature
        print '  metallicity = %.1f' % metallicity
        print '  log gravity = %.1f' % gravity

    return sp

def get_amesdusty_atmosphere(metallicity=0, temperature=5000, gravity=4):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 5000)
    gravity = log gravity (def = 4.0)
    """
    sp = pysynphot.Icat('AMESdusty', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp.flux != 0)[0]
    if len(idx) == 0:
        print 'Could not find AMESdusty Allard+ 2000 atmosphere model for'
        print '  temperature = %d' % temperature
        print '  metallicity = %.1f' % metallicity
        print '  log gravity = %.1f' % gravity

    return sp

def get_phoenix_atmosphere(metallicity=0, temperature=5000, gravity=4):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 5000)
    gravity = log gravity (def = 4.0)
    """
    sp = pysynphot.Icat('phoenix', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp.flux != 0)[0]
    if len(idx) == 0:
        print 'Could not find PHOENIX BT-Settl (Allard+ 2011 atmosphere model for'
        print '  temperature = %d' % temperature
        print '  metallicity = %.1f' % metallicity
        print '  log gravity = %.1f' % gravity

    return sp

def get_cmfgenRot_atmosphere(metallicity=0, temperature=24000, gravity=4.3, rebin=True):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 24000)
    gravity = log gravity (def = 4.3)

    rebin=True: pull from atmospheres at ck04model resolution.
    """
    if rebin:
        sp = pysynphot.Icat('cmfgen_rot_rebin', temperature, metallicity, gravity)
    else:
        sp = pysynphot.Icat('cmfgen_rot', temperature, metallicity, gravity)
        
    # Do some error checking
    idx = np.where(sp.flux != 0)[0]
    if len(idx) == 0:
        print 'Could not find CMFGEN rotating atmosphere model (Fierro+15) for'
        print '  temperature = %d' % temperature
        print '  metallicity = %.1f' % metallicity
        print '  log gravity = %.1f' % gravity

    return sp

def get_cmfgenNoRot_atmosphere(metallicity=0, temperature=22500, gravity=3.98, rebin=True):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 24000)
    gravity = log gravity (def = 4.3)

    rebin=True: pull from atmospheres at ck04model resolution.
    """
    if rebin:
        sp = pysynphot.Icat('cmfgen_norot_rebin', temperature, metallicity, gravity)
    else:
        sp = pysynphot.Icat('cmfgen_norot', temperature, metallicity, gravity)
        
    # Do some error checking
    idx = np.where(sp.flux != 0)[0]
    if len(idx) == 0:
        print 'Could not find CMFGEN rotating atmosphere model (Fierro+15) for'
        print '  temperature = %d' % temperature
        print '  metallicity = %.1f' % metallicity
        print '  log gravity = %.1f' % gravity

    return sp

def get_cmfgenNoRot_atmosphere(metallicity=0, temperature=30000, gravity=4.14):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 30000)
    gravity = log gravity (def = 4.14)
    """
    sp = pysynphot.Icat('cmfgenF15_noRot', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp.flux != 0)[0]
    if len(idx) == 0:
        print 'Could not find CMFGEN non-rotating atmosphere model (Fierro+15) for'
        print '  temperature = %d' % temperature
        print '  metallicity = %.1f' % metallicity
        print '  log gravity = %.1f' % gravity

    return sp

def get_phoenixv16_atmosphere(metallicity=0, temperature=4000, gravity=4, rebin=True):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 5000)
    gravity = log gravity (def = 4.0)

    temp: 2300 - 7000 steps of 100 K; 7000 - 12000 in steps of 200 K
    grav: 0.0 - 6.0, steps of 0.5 (gaurenteed over all temps)

    If rebin = True, pull from spectra that have been rebinned to ck04model resolution;
    this is important for spectrophotometry, otherwise it takes forever
    """
    if (gravity < 0.5):
        print 'Changing gravity for T=', temperature, ' logg=', gravity
        gravity = 0.5
        
    if rebin == True:
        sp = pysynphot.Icat('phoenix_v16_rebin', temperature, metallicity, gravity)
    else:
        sp = pysynphot.Icat('phoenix_v16', temperature, metallicity, gravity)
    
    # Do some error checking
    idx = np.where(sp.flux != 0)[0]
    if len(idx) == 0:
        print 'Could not find PHOENIXv16 (Husser+13) atmosphere model for'
        print '  temperature = %d' % temperature
        print '  metallicity = %.1f' % metallicity
        print '  log gravity = %.1f' % gravity

    return sp

def get_atlas_phoenix_atmosphere(metallicity=0, temperature=4000, gravity=4):
    """
    Return atmosphere that is a linear merge of atlas ck04 model and phoenixV16.

    Only valid for temp of 5250 K, gravity from 0 = 5.0 in steps of 0.5
    """
    sp = pysynphot.Icat('merged', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp.flux != 0)[0]
    if len(idx) == 0:
        print 'Could not find ATLAS-PHOENIX merge atmosphere model for'
        print '  temperature = %d' % temperature
        print '  metallicity = %.1f' % metallicity
        print '  log gravity = %.1f' % gravity

    return sp


#---------OLD MERGED ATMOSPHERES------------#
#def test_merged_atmospheres(metallicity=0, gravity=4):
#    """
#    Compare spectra from Castelli and NextGen at the boundary
#    temperature of 8000 K.
#
#    Compare spectra from NextGen and Phoenix at the boundary
#    temperature of 4000 K.
#    """
#    cast = get_castelli_atmosphere(temperature=8000, 
#                                   metallicity=metallicity, gravity=gravity)
#
#    ngen = get_nextgen_atmosphere(temperature=8000,
#                                  metallicity=metallicity, gravity=gravity)
#
#    # Now Plot the spectra
#    py.figure(1)
#    py.clf()
#    py.loglog(cast.wave, cast.flux, 'r-', label='Castelli')
#    py.plot(ngen.wave, ngen.flux, 'b-', label='NextGen')
#    py.xlabel('Wavelength')
#    py.ylabel('Flux')
#    py.legend()
#    py.xlim(3000, 50000)
#    py.ylim(1e3, 1e8)
#
#
#    ngen = get_nextgen_atmosphere(temperature=4000,
#                                  metallicity=metallicity, gravity=gravity)
#
#    phoe = get_phoenix_atmosphere(temperature=4000,
#                                  metallicity=metallicity, gravity=gravity)
#    # Now Plot the spectra
#    py.figure(2)
#    py.clf()
#    py.loglog(phoe.wave, phoe.flux, 'r-', label='Phoenix')
#    py.plot(ngen.wave, ngen.flux, 'b-', label='NextGen')
#    py.xlabel('Wavelength')
#    py.ylabel('Flux')
#    py.legend()
#    py.xlim(3000, 50000)
#    py.ylim(1, 1e8)
#
#    
#def get_merged_atmosphere(metallicity=0, temperature=20000, gravity=4):
#    if temperature < 4000 or (temperature < 7000 and gravity < 4.0):
#        print 'Phoenix Model Atmosphere Used'
#        return get_phoenix_atmosphere(metallicity=metallicity,
#                                      temperature=temperature,
#                                      gravity=gravity)
#
#    # if temperature < 4000:
#    #     return get_amesdusty_atmosphere(metallicity=metallicity,
#    #                                     temperature=temperature,
#    #                                     gravity=gravity)
#
#    if temperature >= 4000 and temperature < 7000 and gravity >= 4.0:
#        print 'Nextgen atmosphere used'
#        return get_nextgen_atmosphere(metallicity=metallicity,
#                                      temperature=temperature,
#                                      gravity=gravity)
#
#    if temperature >= 7000: 
#        return get_castelli_atmosphere(metallicity=metallicity,
#                                       temperature=temperature,
#                                       gravity=gravity)
#
#---------------------------------------------------------------------#
def get_merged_atmosphere(metallicity=0, temperature=20000, gravity=4):
    """
    If T > 20,000 K : CMFGEN
    20,000 > T > 5500: ATLAS (ck04)
    5500 > T > 5000: ATLAS/PHOENIX merge
    T < 5000: PHEONIXv16 (Husser+13) 
    """
    if temperature < 5000:
        return get_phoenixv16_atmosphere(metallicity=metallicity,
                                      temperature=temperature,
                                      gravity=gravity)

    if (temperature > 5000) & (temperature < 5500):
        #print 'ATLAS and PHOENIX merged atmosphere used'
        #return get_atlas_phoenix_atmosphere(metallicity=metallicity,
        #                                temperature=temperature,
        #                                 gravity=gravity)
        return get_phoenixv16_atmosphere(metallicity=metallicity,
                                      temperature=temperature,
                                      gravity=gravity)
    
    if temperature > 5500:
        return get_castelli_atmosphere(metallicity=metallicity,
                                      temperature=temperature,
                                      gravity=gravity)

    if temperature > 20000:
        return get_castelli_atmosphere(metallicity=metallicity,
                                       temperature=temperature,
                                       gravity=gravity)

        #return get_cmfgenrot_atmosphere(metallicity=metallicity,
        #                               temperature=temperature,
        #                               gravity=gravity)    

#--------------------------------------#
# Atmosphere formatting functions
#--------------------------------------#

def download_CMFGEN_atmospheres(Table_rot, Table_norot):
    """
    Downloads CMFGEN models from
    https://sites.google.com/site/fluxesandcontinuum/home;
    these contain continuum as well as lines.

    Table_rot, Table_norot are tables with the file prefixes
    and model atmosphere parameters, taken by hand from the
    Fierro+15 paper

    Website addresses are hardcoded
    
    Puts downloaded models in the current working directory.
    """
    print 'WARNING: THIS DOES NOT COMPLETELY WORK'
    print '**********************'
    t_rot = Table.read(Table_rot, format='ascii')
    t_norot = Table.read(Table_norot, format='ascii')

    tables = [t_rot, t_norot]
    filenames = [t_rot['col1'], t_norot['col1']]

    # Hardcoded list of webiste addresses
    web_base1 = 'https://sites.google.com/site/fluxesandcontinuum/home/'
    web_base2 = 'https://sites.google.com/site/modelsobmassivestars/'
    web = [web_base1+'009-solar-masses/',web_base1+'012-solar-masses/',
           web_base1+'015-solar-masses/',web_base1+'020-solar-masses/',
           web_base1+'025-solar-masses/',web_base2+'009-solar-masses-tracks/',
           web_base2+'040-solar-masses/',web_base2+'060-solar-masses/',
           web_base1+'085-solar-masses/',web_base1+'120-solar-masses/']
    # Array of masses that matches the website addresses
    mass_arr = np.array([9.,12.,15.,20.,25.,32.,40.,60.,85.,120.])

    # Loop through rotating and unrotating case. First loop is rot, second unrot
    for i in range(2):
        # Extract masses from filenames
        masses = []
        for j in filenames[i]:
            tmp = j.split('m')
            mass = float(tmp[1][:-1])
            masses.append(mass)

        # Download the models webpage by webpage. A bit tricky because masses
        # change slightly within a particular website. THIS IS WHAT FAILS
        for j in range(len(web)):
            if j == 0:
                good = np.where( (masses <= mass_arr[j]) )
            else:
                g = j - 1
                good = np.where( (masses <= mass_arr[j]) &
                                (masses > mass_arr[g]) )
            # Use wget command to pull down the files, and unzip them
            for k in good[0]:
                full = web[j]+'{1:s}.flx.zip'.format(mass_arr[j],filenames[i][k])
                os.system('wget ' + full)
                os.system('unzip '+ filenames[i][k] + '.flx.zip')

    return

def organize_CMFGEN_atmospheres(path_to_dir):
    """
    Organize CMFGEN grid from Fierro+15
    (http://www.astroscu.unam.mx/atlas/index.html)
    into rot and noRot directories

    path_to_dir is from current working directory to directory
    containing the downloaded models. Assumed that models
    and tables describing parameters are in this directory.

    Tables describing parameters MUST be named Table_rot.txt,
    Table_noRot.txt. Made by hand from Tables 3, 4 in Fierro+15.
    These are located in same original directory as atmosphere files

    Will separate files into 2 subdirectories, one rotating and
    the other non-rotating

    *Can't have any other files starting with "t" in model directory to start!*
    """
    # First, record current working directory to return to later
    start_dir = os.getcwd()
    
    # Enter atmosphere directory, collect rotating and non-rotating
    # file names (assumed to all start with "t")
    os.chdir(path_to_dir)
    rot_models = glob.glob("t*r.flx*")
    noRot_models = glob.glob("t*n.flx*")

    # Separate into different subdirectories
    if os.path.exists('cmfgenF15_rot'):
        pass
    else:
        os.mkdir('cmfgenF15_rot')
        os.mkdir('cmfgenF15_noRot')

    for mod in rot_models:
        cmd = 'mv {0:s} cmfgenF15_rot'.format(mod)
        os.system(cmd)

    for mod in noRot_models:
        cmd = 'mv {0:s} cmfgenF15_noRot'.format(mod)
        os.system(cmd)

    # Also move Tables with model parameters into correct directory
    os.system('mv Table_rot.txt cmfgenF15_rot')
    os.system('mv Table_noRot.txt cmfgenF15_noRot')
    
    # Return to original directory
    os.chdir(start_dir)
        
    return

def make_CMFGEN_catalog(path_to_dir):
    """
    Create cdbs catalog.fits of CMFGEN grid from Fierro+15
    (http://www.astroscu.unam.mx/atlas/index.html).
    THIS IS STEP 2, after organize_CMFGEN_atmospheres has
    been run.

    path_to_dir is from current working directory to directory
    containing the rotating or non-rotating models (i.e. cmfgenF15_rot). Also,
    needs to be a Table*.txt file which contains the parameters for all of the
    original models, since params in filename are not precise enough

    Will create catalog.fits file in atmosphere directory with
    description of each model

    *Can't have any other files starting with "t" in model directory to start!*
    """
    # Record current working directory for later
    start_dir = os.getcwd()
    
    # Enter atmosphere directory
    os.chdir(path_to_dir)
   
    # Extract parameters for each atmosphere
    # Note: can't rely on filename for this because not precise enough!!

    #---------OLD: GETTING PARAMS FROM FILENAME-------#
    # Collect file names (assumed to all start with "t")
    #files = glob.glob("t*")
    #for name in files:
    #    tmp = name.split('l')
    #    temp = float(tmp[0][1:]) * 100.0 # In kelvin

    #    lumtmp = tmp[1].split('_')
    #    lum = float(lumtmp[0][:-5]) * 1000.0 # In L_sun

    #    mass = float(lumtmp[0][5:-1]) # In M_sun
        
        # Need to calculate log g from T and L (cgs)
    #    lum_sun = 3.846 * 10**33 # erg/s
    #    M_sun = 2 * 10**33 # g
    #    G_si = 6.67 * 10**(-8) # cgs
    #    sigma_si = 5.67 * 10**(-5) # cgs

    #    g = (G_si * mass * M_sun * 4 * np.pi * sigma_si * temp**4) / \
    #      (lum * lum_sun)
    #    logg = np.log10(g)
    #---------------------------------------------------#

    # Read table with atmosphere params
    table = glob.glob('Table_*')
    t = Table.read(table[0], format = 'ascii')
    names = t['col1']
    temps = t['col2']
    logg = t['col4']

    # Create catalog.fits file
    index_str = []
    name_str = []
    for i in range(len(names)):
        index = '{0:5.0f},0.0,{1:3.2f}'.format(temps[i], logg[i])

        #---NOTE: THE FOLLOWING DEPENDS ON FINAL LOCATION OF CATALOG FILE---#
        #path = path_to_dir + '/' + names[i]
        path = names[i] + '.fits[Flux]'
        
        index_str.append(index)
        name_str.append(path)
    
    catalog = Table([index_str, name_str], names = ('INDEX', 'FILENAME'))

    # Create catalog.fits file in directory with the models
    catalog.write('catalog.fits', format = 'fits')
    
    # Move back to original directory, create the catalog.fits file
    os.chdir(start_dir)
    
    return

def cdbs_cmfgen(path_to_dir, path_to_cdbs_dir):
    """
    Code to put cmfgen models into cdbs format and adds proper unit keyword in
    fits header. Save as fits file

    path_to_dir goes from current directory to cmfgen_rot or cmfgen_norot
    directory with the *.flx models. Note that these files have already been
    organized using organize_CMFGEN_atmospheres code.

    path_to_cdbs_dir goes from current directory to cdbs/grid/cmfgen_rot or
    cmfgen_norot directory. Will copy new fits files to this directory.
    This directory must already exist!
    """
    # Save starting directory for later, move into path_to_dir directory
    start_dir = os.getcwd()
    os.chdir(path_to_dir)

    # Collect the filenames, make necessary changes to each one
    files = glob.glob('*.flx')

    # Need to make brand-new fits tables with data we want.
    counter = 0
    for i in files:
        counter += 1
        # Open file, extract useful info
        t = Table.read(i, format='ascii')
        wave = t['col1']
        flux = t['col2'] # Flux is already in erg/cm^2/s/A

        # Need to eliminate duplicate entries (pysynphot crashes)
        unique = np.unique(wave, return_index=True)
        wave = wave[unique[1]]
        flux = flux[unique[1]]
        
        # Make fits table from individual columns. 
        c0 = fits.Column(name='Wavelength', format='D', array=wave)
        c1 = fits.Column(name='Flux', format='E', array=flux)

        cols = fits.ColDefs([c0, c1])
        tbhdu = fits.BinTableHDU.from_columns(cols)

        #Adding unit keywords        
        tbhdu.header['TUNIT1'] = 'ANGSTROM'
        tbhdu.header['TUNIT2'] = 'FLAM'

        prihdu = fits.PrimaryHDU()
    
        finalhdu = fits.HDUList([prihdu, tbhdu])
        finalhdu.writeto(i[:-4]+'.fits', clobber=True)
        
        print 'Done {0:2.0f} of {1:2.0f}'.format(counter, len(files))

    # Return to original directory, copy over new .fits files to cdbs directory
    os.chdir(start_dir)
    cmd = 'mv {0:s}/*.fits {1:s}'.format(path_to_dir, path_to_cdbs_dir)
    os.system(cmd)

    return

def rebin_cmfgen(cdbs_path, rot=True):
    """
    Rebin cmfgen_rot and cmfgen_norot models to atlas ck04 resolution;
    this makes spectrophotometry MUCH faster

    cdbs_path: path to cdbs directory
    rot=True for rotating models (cmfgen_rot), False for non-rotating models
    
    makes new directory in cdbs/grid: cmfgen_rot_rebin or cmfgen_norot_rebin
    """
    # Get an atlas ck04 model, we will use this to set wavelength grid
    sp_atlas = get_castelli_atmosphere()

    # Open a fits table for an existing cmfgen model; we will steal the header.
    # Also define paths to new rebin directories
    if rot == True:
        tmp = cdbs_path+'/grid/cmfgen_rot/t0200l0008m009r.fits'
        path = cdbs_path+'/grid/cmfgen_rot_rebin/'
        orig_path = cdbs_path+'/grid/cmfgen_rot/'
    else:
        tmp = cdbs_path+'/grid/cmfgen_norot/t0200l0007m009n.fits'
        path = cdbs_path+'/grid/cmfgen_norot_rebin/'
        orig_path = cdbs_path+'/grid/cmfgen_norot/'
        
    cmfgen_hdu = fits.open(tmp)
    header0 = cmfgen_hdu[0].header
    # Create rebin directories if they don't already exist. Copy over
    # catalog.fits file from original directory (will be the same)
    if not os.path.exists(path):
        os.mkdir(path)
        cmd = 'cp {0:s}catalog.fits {1:s}'.format(orig_path, path)
        os.system(cmd)

    # Read in the catalog.fits file
    cat = fits.getdata(orig_path + 'catalog.fits')
    files_all = [cat[ii][1].split('[')[0] for ii in range(len(cat))]

    # First column in new files will be for [atlas] wavelength
    c0 = fits.Column(name='Wavelength', format='D', array=sp_atlas.wave)   

    # For each catalog.fits entry, read the unbinned spectrum and rebin to
    # the atlas resolution. Make a new fits file in rebin directory
    count = 0
    for ff in range(len(files_all)):
        count += 1
        # Extract the temp, Z, logg
        vals = cat[ff][0].split(',')
        temp = float(vals[0])
        metal = float(vals[1])
        grav = float(vals[2])
    
        # Fetch the spectrum
        if rot == True:           
            sp = pysynphot.Icat('cmfgen_rot', temp, metal, grav)
        else:
            sp = pysynphot.Icat('cmfgen_norot', temp, metal, grav)

        # Rebin
        flux_rebin = rebin_spec(sp.wave, sp.flux, sp_atlas.wave)
        c1 = fits.Column(name='Flux', format='E', array=flux_rebin)                

        # Make the FITS file from the columns with header
        cols = fits.ColDefs([c0,c1])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        prihdu = fits.PrimaryHDU(header=header0)
        tbhdu.header['TUNIT1'] = 'ANGSTROM'
        tbhdu.header['TUNIT2'] = 'FLAM'

        # Write hdu to new directory with same filename
        finalhdu = fits.HDUList([prihdu, tbhdu])
        finalhdu.writeto(path+files_all[ff])

        print 'Finished file {0} of {1}'.format(count, len(files_all))            
    return


def organize_PHOENIXv16_atmospheres(path_to_dir):
    """
    Construct the Phoenix Husser+13 atmopsheres for each model. Combines the
    fluxes from the *HiRES.fits files and the wavelengths of the
    WAVE_PHONEIX-ACES-AGSS-COND-2011.fits file.

    path_to_dir is the path to the directory containing all of the downloaded
    files

    Creates new fits files for each atmosphere: phoenix_mm00_<temp>.fits,
    which contains columns for the log g (column header = g#.#). Puts
    atmospheres in new directory phoenixm00
    """
    # Save current directory for return later, move into working dir
    start_dir = os.getcwd()
    os.chdir(path_to_dir)

    # If it doesn't already exist, create the phoenixm00 subdirectory
    if os.path.exists('phoenixm00'):
        pass
    else:
        os.mkdir('phoenixm00')
    
    # Extract wavelength array, make column for later
    wavefile = fits.open('WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
    wave = wavefile[0].data
    wavefile.close()
    wave_col = Column(wave, name = 'WAVELENGTH')

    # Create temp array for Husser+13 grid (given in paper)
    temp_arr = np.arange(2300, 7001, 100)
    temp_arr = np.append(temp_arr, np.arange(7000, 12001, 200))

    print 'Looping though all temps'
    # For each temp, build file containing the flux for all gravities
    i = 0
    for temp in temp_arr:
        files = glob.glob('lte{0:05d}-*-HiRes.fits'.format(temp))
        # Start the table with the wavelength column
        t = Table()
        t.add_column(wave_col)
        for f in files:
            # Extract the logg out of filename
            tmp = f.split('-')
            logg = tmp[1]
            
            # Extract fluxes from file
            spectrum = fits.open(f)
            flux = spectrum[0].data
            spectrum.close()

            # Make Column object with fluxes, add to table
            col = Column(flux, name = 'g{0:2.1f}'.format(float(logg)))
            t.add_column(col)
            
        # Now, construct final fits file for the given temp
        outname = 'phoenixm00_{0:05d}.fits'.format(temp)
        t.write('phoenixm00/' + outname, format = 'fits', overwrite = True) 
        
        # Progress counter for user
        i += 1
        print 'Done {0:d} of {1:d}'.format(i, len(temp_arr))

    # Return to original directory
    os.chdir(start_dir)
    return

def make_PHOENIXv16_catalog(path_to_dir):
    """
    Makes catalog.fits file for Husser+13 phoenix models. Assumes that
    construct_atmospheres has been run already, and that the models lie
    in subdirectory phoenixm00.

    path_to_directory is the path to the directory with the reformatted
    models (i.e. the output from construct_atmospheres, phoenixm00)
    
    Puts catalog.fits file in directory the user starts in
    """
    # Save starting directory for later, move into working directory
    start_dir = os.getcwd()
    os.chdir(path_to_dir)

    # Collect the filenames. Each is a unique temp with many different log g's
    files = glob.glob('*.fits')

    # Create the catalog.fits file, row by row
    index_arr = []
    filename_arr = []
    for i in files:
        # Get log g values from the column header the file
        t = Table.read(i, format='fits')
        keys = t.keys()
        logg_vals = keys[1:]
        
        # Extract temp from filename
        name = i.split('_')
        temp = float(name[1][:-5])
        for j in logg_vals:
            logg = float(j[1:])
            index = '{0:5.0f},0.0,{1:2.1f}'.format(temp, logg)
            filename = path_to_dir + i + '[' + j + ']'
            # Add row to table
            index_arr.append(index)
            filename_arr.append(filename)

    catalog = Table([index_arr, filename_arr], names=('INDEX', 'FILENAME'))
    
    # Return to starting directory, write catalog
    os.chdir(start_dir)
    catalog.write('catalog.fits', format='fits', overwrite=True)
    
    return

def cdbs_PHOENIXv16(path_to_cdbs_dir):
    """
    Put the PHOENIXv16 (Husser+13) fits files into cdbs format. This primarily
    consists of adjusting the flux units from [erg/s/cm^2/cm] to [erg/s/cm^2/A]
    and adding the appropriate keywords to the fits header.

    path_to_cdbs_dir goes from current working directory to phoenixm00 directory
    in cdbs/grids/phoenix_v16. Note that these files have already been organized
    using orgnaize_PHOENIXv16_atmospheres code.

    Overwrites original files in directory
    """
    # Save starting directory for later, move into working directory
    start_dir = os.getcwd()
    os.chdir(path_to_cdbs_dir)

    # Collect the filenames, make necessary changes to each one
    files = glob.glob('*.fits')

    # Need to make brand-new fits tables with data we want.
    counter = 0
    for i in files:
        counter += 1
        # Open file, extract useful info
        hdu = fits.open(i)
        header_0 = hdu[0].header
        header_1 = hdu[1].header
        sci = hdu[1].data

        # Remake fits table from individual columns, multiplying each flux
        # column by 10^-8 for conversion
        # This gets messy due to changing number of columns
        c0 = fits.Column(name='Wavelength', format='D', array=sci.field(0))
        # This particular column only exists for lower temp models
        if counter <= 34:
            c1 = fits.Column(name='g0.0', format='E', array=sci.field(1)*10**-8)
            c2 = fits.Column(name='g0.5', format='E', array=sci.field(2)*10**-8)
            c3 = fits.Column(name='g1.0', format='E', array=sci.field(3)*10**-8)
            c4 = fits.Column(name='g1.5', format='E', array=sci.field(4)*10**-8)
            c5 = fits.Column(name='g2.0', format='E', array=sci.field(5)*10**-8)
            c6 = fits.Column(name='g2.5', format='E', array=sci.field(6)*10**-8)
            c7 = fits.Column(name='g3.0', format='E', array=sci.field(7)*10**-8)
            c8 = fits.Column(name='g3.5', format='E', array=sci.field(8)*10**-8)
            c9 = fits.Column(name='g4.0', format='E', array=sci.field(9)*10**-8)
            c10 = fits.Column(name='g4.5', format='E', array=sci.field(10)*10**-8)
            c11 = fits.Column(name='g5.0', format='E', array=sci.field(11)*10**-8)
            c12 = fits.Column(name='g5.5', format='E', array=sci.field(12)*10**-8)
            c13 = fits.Column(name='g6.0', format='E', array=sci.field(13)*10**-8)
        elif counter <= 37:
            c2 = fits.Column(name='g0.5', format='E', array=sci.field(1)*10**-8)
            c3 = fits.Column(name='g1.0', format='E', array=sci.field(2)*10**-8)
            c4 = fits.Column(name='g1.5', format='E', array=sci.field(3)*10**-8)
            c5 = fits.Column(name='g2.0', format='E', array=sci.field(4)*10**-8)
            c6 = fits.Column(name='g2.5', format='E', array=sci.field(5)*10**-8)
            c7 = fits.Column(name='g3.0', format='E', array=sci.field(6)*10**-8)
            c8 = fits.Column(name='g3.5', format='E', array=sci.field(7)*10**-8)
            c9 = fits.Column(name='g4.0', format='E', array=sci.field(8)*10**-8)
            c10 = fits.Column(name='g4.5', format='E', array=sci.field(9)*10**-8)
            c11 = fits.Column(name='g5.0', format='E', array=sci.field(10)*10**-8)
            c12 = fits.Column(name='g5.5', format='E', array=sci.field(11)*10**-8)
            c13 = fits.Column(name='g6.0', format='E', array=sci.field(12)*10**-8)
        elif counter <= 54:
            c3 = fits.Column(name='g1.0', format='E', array=sci.field(1)*10**-8)
            c4 = fits.Column(name='g1.5', format='E', array=sci.field(2)*10**-8)
            c5 = fits.Column(name='g2.0', format='E', array=sci.field(3)*10**-8)
            c6 = fits.Column(name='g2.5', format='E', array=sci.field(4)*10**-8)
            c7 = fits.Column(name='g3.0', format='E', array=sci.field(5)*10**-8)
            c8 = fits.Column(name='g3.5', format='E', array=sci.field(6)*10**-8)
            c9 = fits.Column(name='g4.0', format='E', array=sci.field(7)*10**-8)
            c10 = fits.Column(name='g4.5', format='E', array=sci.field(8)*10**-8)
            c11 = fits.Column(name='g5.0', format='E', array=sci.field(9)*10**-8)
            c12 = fits.Column(name='g5.5', format='E', array=sci.field(10)*10**-8)
            c13 = fits.Column(name='g6.0', format='E', array=sci.field(11)*10**-8)
        elif counter <= 59:
            c4 = fits.Column(name='g1.5', format='E', array=sci.field(1)*10**-8)
            c5 = fits.Column(name='g2.0', format='E', array=sci.field(2)*10**-8)
            c6 = fits.Column(name='g2.5', format='E', array=sci.field(3)*10**-8)
            c7 = fits.Column(name='g3.0', format='E', array=sci.field(4)*10**-8)
            c8 = fits.Column(name='g3.5', format='E', array=sci.field(5)*10**-8)
            c9 = fits.Column(name='g4.0', format='E', array=sci.field(6)*10**-8)
            c10 = fits.Column(name='g4.5', format='E', array=sci.field(7)*10**-8)
            c11 = fits.Column(name='g5.0', format='E', array=sci.field(8)*10**-8)
            c12 = fits.Column(name='g5.5', format='E', array=sci.field(9)*10**-8)
            c13 = fits.Column(name='g6.0', format='E', array=sci.field(10)*10**-8)
        else:
            c5 = fits.Column(name='g2.0', format='E', array=sci.field(1)*10**-8)
            c6 = fits.Column(name='g2.5', format='E', array=sci.field(2)*10**-8)
            c7 = fits.Column(name='g3.0', format='E', array=sci.field(3)*10**-8)
            c8 = fits.Column(name='g3.5', format='E', array=sci.field(4)*10**-8)
            c9 = fits.Column(name='g4.0', format='E', array=sci.field(5)*10**-8)
            c10 = fits.Column(name='g4.5', format='E', array=sci.field(6)*10**-8)
            c11 = fits.Column(name='g5.0', format='E', array=sci.field(7)*10**-8)
            c12 = fits.Column(name='g5.5', format='E', array=sci.field(8)*10**-8)
            c13 = fits.Column(name='g6.0', format='E', array=sci.field(9)*10**-8)
                        
        if counter <= 35:
            cols = fits.ColDefs([c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13])
        elif counter <= 37:
            cols = fits.ColDefs([c0,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13])
        elif counter <= 54:
            cols = fits.ColDefs([c0,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13])
        elif counter <= 59:
            cols = fits.ColDefs([c0,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13])
        else:
            cols = fits.ColDefs([c0,c5,c6,c7,c8,c9,c10,c11,c12,c13])
        tbhdu = fits.BinTableHDU.from_columns(cols)

        # Copying over the older headers, adding unit keywords
        prihdu = fits.PrimaryHDU(header=header_0)
        tbhdu.header['TUNIT1'] = 'ANGSTROM'
        tbhdu.header['TUNIT2'] = 'FLAM'
        tbhdu.header['TUNIT3'] = 'FLAM'
        tbhdu.header['TUNIT4'] = 'FLAM'
        tbhdu.header['TUNIT5'] = 'FLAM'
        tbhdu.header['TUNIT6'] = 'FLAM'
        tbhdu.header['TUNIT7'] = 'FLAM'
        tbhdu.header['TUNIT8'] = 'FLAM'
        tbhdu.header['TUNIT9'] = 'FLAM'
        tbhdu.header['TUNIT10'] = 'FLAM'
        tbhdu.header['TUNIT11'] = 'FLAM'
        tbhdu.header['TUNIT12'] = 'FLAM'
        tbhdu.header['TUNIT13'] = 'FLAM'
        tbhdu.header['TUNIT14'] = 'FLAM'
    
        finalhdu = fits.HDUList([prihdu, tbhdu])
        finalhdu.writeto(i, clobber=True)

        hdu.close()
        print 'Done {0:2.0f} of {1:2.0f}'.format(counter, len(files))

    return

def rebin_phoenixV16(cdbs_path):
    """
    Rebin phoenixV16 models to atlas ck04 resolution; this makes
    spectrophotometry MUCH faster

    makes new directory in cdbs/grid: phoenix_v16_rebin

    cdbs_path: path to cdbs directory
    """
    # Get an atlas ck04 model, we will use this to set wavelength grid
    sp_atlas = get_castelli_atmosphere()

    # Open a fits table for an existing phoenix model; we will steal the header
    tmp = cdbs_path+'/grid/phoenix_v16/phoenixm00/phoenixm00_02400.fits'
    phoenix_hdu = fits.open(tmp)
    header0 = phoenix_hdu[0].header

    # Create cdbs/grid directory for rebinned models
    path = cdbs_path+'/grid/phoenix_v16_rebin/'
    if not os.path.exists(path):
        os.mkdir(path)
        os.mkdir(path+'phoenixm00')

    # Read in the existing catalog.fits file and rebin every spectrum.
    cat = fits.getdata(cdbs_path + '/grid/phoenix_v16/catalog.fits')
    files_all = [cat[ii][1].split('[')[0] for ii in range(len(cat))]
    temp_arr = np.zeros(len(files_all), dtype=float)
    logg_arr = np.zeros(len(files_all), dtype=float)
    metal_arr = np.zeros(len(files_all), dtype=float)

    for ff in range(len(files_all)):
        vals = cat[ff][0].split(',')

        temp_arr[ff] = float(vals[0])
        metal_arr[ff] = float(vals[1])
        logg_arr[ff] = float(vals[2])

    metal_uniq = np.unique(metal_arr)
    temp_uniq = np.unique(temp_arr)

    for mm in range(len(metal_uniq)):
        metal = metal_uniq[mm] # metallicity
        
        for tt in range(len(temp_uniq)):
            temp = temp_uniq[tt] # temperature

            # Pick out the list of gravities for this T, Z combo            
            idx = np.where((metal_arr == metal) & (temp_arr == temp))[0]
            logg_exist = logg_arr[idx]

            # All gravities will go in one file. Here is the output
            # file name.
            outfile = path + files_all[idx[0]].split('[')[0]
            
            # Build a columns array. One column for each gravity.
            cols_arr = []

            # Make the wavelength column, which is first in the cols array.
            c0 = fits.Column(name='Wavelength', format='D', array=sp_atlas.wave)
            cols_arr.append(c0)
            
            for gg in range(len(logg_exist)):
                grav = logg_exist[gg] # gravity

                # Fetch the spectrum                
                sp = pysynphot.Icat('phoenix_v16', temp, metal, grav)
                flux_rebin = rebin_spec(sp.wave, sp.flux, sp_atlas.wave)

                # Store the spectrum
                name = 'g{0:3.1f}'.format(grav)
                col = fits.Column(name=name, format='E', array=flux_rebin)
                cols_arr.append(col)
                

            # Make the FITS file from the columns with header.
            cols = fits.ColDefs(cols_arr)
            tbhdu = fits.BinTableHDU.from_columns(cols)
            prihdu = fits.PrimaryHDU(header=header0)
            tbhdu.header['TUNIT1'] = 'ANGSTROM'
            for gg in range(len(logg_exist)):
                tbhdu.header['TUNIT{0:d}'.format(gg+2)] = 'FLAM'

            # Write hdu
            finalhdu = fits.HDUList([prihdu, tbhdu])
            # don't have clobber to protect original files.
            finalhdu.writeto(path + outfile)   

            print 'Finished file ' + outfile + ' with gravities: ', logg_exist
            

    return


def rebin_spec(wave, specin, wavnew):
    """
    Helper routine to rebin spectra. TAKEN FROM ASTROBETTER BLOG FROM JESSICA:
    http://www.astrobetter.com/blog/2013/08/12/
    python-tip-re-sampling-spectra-with-pysynphot/
    """
    spec = pysynphot.spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = pysynphot.spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = pysynphot.observation.Observation(spec, filt, binset=wavnew, force='taper')
 
    return obs.binflux
