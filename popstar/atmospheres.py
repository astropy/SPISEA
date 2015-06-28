import logging
import numpy as np
import pysynphot

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

                  2600 -  4000            100 
                  4000 - 10000            200 

    log gravity (def = 4.0) in the range of 3.5 - 6.0 in 0.5 increments
    """
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

def test_merged_atmospheres(metallicity=0, gravity=4):
    """
    Compare spectra from Castelli and NextGen at the boundary
    temperature of 8000 K.

    Compare spectra from NextGen and Phoenix at the boundary
    temperature of 4000 K.
    """
    cast = get_castelli_atmosphere(temperature=8000, 
                                   metallicity=metallicity, gravity=gravity)

    ngen = get_nextgen_atmosphere(temperature=8000,
                                  metallicity=metallicity, gravity=gravity)

    # Now Plot the spectra
    py.figure(1)
    py.clf()
    py.loglog(cast.wave, cast.flux, 'r-', label='Castelli')
    py.plot(ngen.wave, ngen.flux, 'b-', label='NextGen')
    py.xlabel('Wavelength')
    py.ylabel('Flux')
    py.legend()
    py.xlim(3000, 50000)
    py.ylim(1e3, 1e8)


    ngen = get_nextgen_atmosphere(temperature=4000,
                                  metallicity=metallicity, gravity=gravity)

    phoe = get_phoenix_atmosphere(temperature=4000,
                                  metallicity=metallicity, gravity=gravity)
    # Now Plot the spectra
    py.figure(2)
    py.clf()
    py.loglog(phoe.wave, phoe.flux, 'r-', label='Phoenix')
    py.plot(ngen.wave, ngen.flux, 'b-', label='NextGen')
    py.xlabel('Wavelength')
    py.ylabel('Flux')
    py.legend()
    py.xlim(3000, 50000)
    py.ylim(1, 1e8)

    


def get_merged_atmosphere(metallicity=0, temperature=20000, gravity=4):
    if temperature < 4000 or (temperature < 7000 and gravity < 4.0):
        print 'Phoenix Model Atmosphere Used'
        return get_phoenix_atmosphere(metallicity=metallicity,
                                      temperature=temperature,
                                      gravity=gravity)

    # if temperature < 4000:
    #     return get_amesdusty_atmosphere(metallicity=metallicity,
    #                                     temperature=temperature,
    #                                     gravity=gravity)

    if temperature >= 4000 and temperature < 7000 and gravity >= 4.0:
        print 'Nextgen atmosphere used'
        return get_nextgen_atmosphere(metallicity=metallicity,
                                      temperature=temperature,
                                      gravity=gravity)

    if temperature >= 7000: 
        return get_castelli_atmosphere(metallicity=metallicity,
                                       temperature=temperature,
                                       gravity=gravity)


        


