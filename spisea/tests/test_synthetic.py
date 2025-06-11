import time
import numpy as np
import pylab as plt
import numpy as np
from spisea import reddening, evolution, atmospheres, ifmr
from spisea import synthetic as syn
from spisea.imf import imf
from spisea.imf import multiplicity
import pysynphot
import os
import pdb
from scipy.spatial import cKDTree as KDTree

def test_isochrone(plot=False):
    logAge = 6.7
    AKs = 2.7
    distance = 4000

    startTime = time.time()
    iso = syn.Isochrone(logAge, AKs, distance)
    print('Test completed in: %d seconds' % (time.time() - startTime))
    # Typically takes 104 - 120 seconds.
    # Limited by pysynphot.Icat call in atmospheres.py

    assert iso.points.meta['LOGAGE'] == logAge
    assert iso.points.meta['AKS'] == AKs
    assert iso.points.meta['DISTANCE'] == distance
    assert len(iso.points) > 100
    
    if plot:
        plt.figure(1) 
        iso.plot_HR_diagram()
        
        plt.figure(2)
        iso.plot_mass_luminosity()

    return iso

def test_iso_wave():
    """
    Test to make sure isochrones generated have spectra with the proper 
    wavelength range, and that the user has control over that wavelength
    range (propagated through IsochronePhot)
    """
    # Define isochrone parameters
    logAge = np.log10(5*10**6.) # Age in log(years)
    AKs = 0.8 # extinction in mags
    dist = 4000 # distance in parsec

    # Define evolution/atmosphere models and extinction law (optional)
    evo_model = evolution.MergedBaraffePisaEkstromParsec() 
    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawHosek18b()

    # Also specify filters for synthetic photometry (optional). Here we use 
    # the HST WFC3-IR F127M, F139M, and F153M filters
    filt_list = ['wfc3,ir,f127m']

    # First, let's make sure the vega spectrum has the proper limits
    vega = syn.Vega()

    assert np.min(vega.wave) == 995
    assert np.max(vega.wave) == 100200

    # Make Isochrone object. Will use wave_range = [3000,52000].
    # Make sure range matches to resolution of atmosphere.
    wave_range1 = [3000, 52000]
    my_iso = syn.IsochronePhot(logAge, AKs, dist,
                            evo_model=evo_model, atm_func=atm_func,
                            red_law=red_law, filters=filt_list,
                            mass_sampling=10, wave_range=wave_range1,
                            recomp=True)

    test = my_iso.spec_list[0]

    assert np.min(test.wave) == 3010
    assert np.max(test.wave) == 51900

    # Now let's try changing the wave range. Is it carried through
    # properly?
    wave_range2 = [1200, 90000]
    my_iso = syn.IsochronePhot(logAge, AKs, dist,
                            evo_model=evo_model, atm_func=atm_func,
                            red_law=red_law, filters=filt_list,
                            mass_sampling=10, wave_range=wave_range2,
                            recomp=True)

    test2 = my_iso.spec_list[0]

    assert np.min(test2.wave) == 1205
    assert np.max(test2.wave) == 89800

    # Does the error exception catch the bad wave_range?
    wave_range3 = [1200, 1000000]
    try:
        my_iso = syn.IsochronePhot(logAge, AKs, dist,
                                evo_model=evo_model, atm_func=atm_func,
                                red_law=red_law, filters=filt_list,
                                mass_sampling=10, wave_range=wave_range3,
                                recomp=True)
        print('WAVE TEST FAILED!!! Should have crashed here, wavelength range out of bounds')
        raise ValueError() 
    except:
        print('Wavelength out of bound condition passed. Test is good')
        pass
    return

def test_IsochronePhot(plot=False):
    logAge = 6.7
    AKs = 2.7
    distance = 4000
    filt_list = ['wfc3,ir,f127m', 'nirc2,J']
    mass_sampling=1
    iso_dir = 'iso/'

    evo_model = evolution.MISTv1()
    atm_func = atmospheres.get_merged_atmosphere
    redlaw = reddening.RedLawNishiyama09()

    startTime = time.time()
    iso = syn.IsochronePhot(logAge, AKs, distance, evo_model=evo_model,
                                atm_func=atm_func, red_law=redlaw,
                                filters=filt_list,
                                mass_sampling=mass_sampling, iso_dir=iso_dir)
    endTime = time.time()
    print('IsochronePhot generated in: %d seconds' % (endTime - startTime))
    # Typically takes 120 seconds if file is regenerated.
    # Limited by pysynphot.Icat call in atmospheres.py

    assert iso.points.meta['LOGAGE'] == logAge
    assert iso.points.meta['AKS'] == AKs
    assert iso.points.meta['DISTANCE'] == distance
    assert len(iso.points) > 100

    assert 'm_nirc2_J' in iso.points.colnames

    if plot:
        plt.figure(1) 
        iso.plot_CMD('mag814w', 'mag160w')
        
        plt.figure(2)
        iso.plot_mass_magnitude('mag160w')

    # Finally, let's test the isochronePhot file generation
    assert os.path.exists('{0}/iso_{1:.2f}_{2:4.2f}_{3:4s}_p00.fits'.format(iso_dir, logAge,
                                                                                AKs, str(distance).zfill(5)))
    
    # Check 1: If we try to remake the isochrone, does it read the file rather than
    # making a new one
    iso_new = syn.IsochronePhot(logAge, AKs, distance, evo_model=evo_model,
                                atm_func=atm_func, red_law=redlaw,
                                filters=filt_list,
                                mass_sampling=mass_sampling, iso_dir=iso_dir)

    assert iso_new.recalc == False
    
    # Check 2: Confirm that adding a new column to an existing isochrone works properly.
    #    Does the new filter get added to the isochrone? And the old ones still there?
    #    Does the computed data for the new filter match the same result if you fully regenerate the isochrone?
    iso_new_addfilt = syn.IsochronePhot(logAge, AKs, distance, evo_model=evo_model,
                                atm_func=atm_func, red_law=redlaw,
                                filters=filt_list+['2mass,Ks'],
                                mass_sampling=mass_sampling, iso_dir=iso_dir)

    assert iso_new_addfilt.recalc == False
    assert 'm_2mass_Ks' in iso_new_addfilt.points.colnames
    assert 'm_nirc2_J' in iso_new_addfilt.points.colnames
    
    iso_new_3filt = syn.IsochronePhot(logAge, AKs, distance, evo_model=evo_model,
                                atm_func=atm_func, red_law=redlaw,
                                filters=filt_list+['2mass,Ks'],
                                mass_sampling=mass_sampling, iso_dir=iso_dir,
                                recomp=True)
    np.testing.assert_almost_equal(iso_new_addfilt.points['m_2mass_Ks'], iso_new_3filt.points['m_2mass_Ks'])
    assert iso_new_3filt.recalc==True

    # Check 3: If we change evo model, atmo model, or redlaw,
    # does IsochronePhot regenerate the isochrone and overwrite the existing one?
    evo2 = evolution.MergedBaraffePisaEkstromParsec()
    mass_sampling=20

    iso_new = syn.IsochronePhot(logAge, AKs, distance, evo_model=evo2,
                                atm_func=atm_func, red_law=redlaw,
                                filters=filt_list,
                                mass_sampling=mass_sampling, iso_dir=iso_dir)

    assert iso_new.recalc == True

    redlaw2 = reddening.RedLawHosek18b()
    iso_new = syn.IsochronePhot(logAge, AKs, distance, evo_model=evo2,
                                atm_func=atm_func, red_law=redlaw2,
                                filters=filt_list,
                                mass_sampling=mass_sampling, iso_dir=iso_dir)

    assert iso_new.recalc == True

    atm2 = atmospheres.get_castelli_atmosphere
    iso_new = syn.IsochronePhot(logAge, AKs, distance, evo_model=evo2,
                                atm_func=atm2, red_law=redlaw2,
                                filters=filt_list,
                                mass_sampling=mass_sampling, iso_dir=iso_dir)

    assert iso_new.recalc == True

    return

def test_ResolvedCluster():
    # Define cluster parameters
    logAge = 6.7
    AKs = 2.4
    distance = 4000
    cluster_mass = 10**5.
    mass_sampling=5

    # Test filters
    filt_list = ['nirc2,J', 'nirc2,Kp']

    startTime = time.time()
    
    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere

    red_law = reddening.RedLawNishiyama09()
    
    iso = syn.IsochronePhot(logAge, AKs, distance,
                            evo_model=evo, atm_func=atm_func,
                            red_law=red_law, filters=filt_list,
                            mass_sampling=mass_sampling)

    print('Constructed isochrone: %d seconds' % (time.time() - startTime))

    # Now to create the cluster.
    imf_mass_limits = np.array([0.07, 0.5, 1, np.inf])
    imf_powers = np.array([-1.3, -2.3, -2.3])

    ##########
    # Start without multiplicity
    ##########
    my_imf1 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=None)
    print('Constructed IMF: %d seconds' % (time.time() - startTime))
    
    cluster1 = syn.ResolvedCluster(iso, my_imf1, cluster_mass)
    clust1 = cluster1.star_systems
    print('Constructed cluster: %d seconds' % (time.time() - startTime))

    # Check that stars are returned
    assert len(clust1) > 0

    # Check that the total mass in stars is less than requested (no compact objects).
    cluster_mass_out = clust1['systemMass'].sum()
    assert cluster_mass_out < cluster_mass

    plt.figure(3)
    plt.clf()
    plt.plot(clust1['m_nirc2_J'] - clust1['m_nirc2_Kp'], clust1['m_nirc2_J'], 'r.')
    plt.plot(iso.points['m_nirc2_J'] - iso.points['m_nirc2_Kp'], iso.points['m_nirc2_J'], 'c.')
    plt.gca().invert_yaxis()

    # *** Visual Inspections: ***
    #  - check that points (red) fall between isochrone points (blue)

    ##########
    # Test with multiplicity
    ##########
    multi = multiplicity.MultiplicityUnresolved()
    my_imf2 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=multi)
    print('Constructed IMF with multiples: %d seconds' % (time.time() - startTime))
    
    cluster2 = syn.ResolvedCluster(iso, my_imf2, cluster_mass)
    clust2 = cluster2.star_systems
    print('Constructed cluster with multiples: %d seconds' % (time.time() - startTime))

    assert len(clust2) > 0
    assert len(cluster2.companions) > 0
    assert np.sum(clust2['N_companions']) == len(cluster2.companions)

    ##########
    # Plots 
    ##########
    # Plot an IR CMD and compare cluster members to isochrone.
    plt.figure(1)
    plt.clf()
    plt.plot(clust1['m_nirc2_J'] - clust1['m_nirc2_Kp'], clust1['m_nirc2_J'], 'r.')
    plt.plot(clust2['m_nirc2_J'] - clust2['m_nirc2_Kp'], clust2['m_nirc2_J'], 'b.')
    plt.plot(iso.points['m_nirc2_J'] - iso.points['m_nirc2_Kp'], iso.points['m_nirc2_J'], 'c-')
    plt.gca().invert_yaxis()
    plt.xlabel('J - Kp (mag)')
    plt.ylabel('J (mag')

    # Plot a mass-magnitude relationship.
    plt.figure(2)
    plt.clf()
    plt.semilogx(clust1['mass'], clust1['m_nirc2_J'], 'r.')
    plt.semilogx(clust2['mass'], clust2['m_nirc2_J'], 'r.')
    plt.gca().invert_yaxis()
    plt.xlabel('Mass (Msun)')
    plt.ylabel('J (mag)')
    
    # # Plot the spectrum of the most massive star
    # idx = cluster.mass.argmax()
    # plt.clf()
    # plt.plot(cluster.stars[idx].wave, cluster.stars[idx].flux, 'k.')

    # # Plot an integrated spectrum of the whole cluster.
    # wave, flux = cluster.get_integrated_spectrum()
    # plt.clf()
    # plt.plot(wave, flux, 'k.')

    return

def test_ResolvedClusterDiffRedden():
    logAge = 6.7
    AKs = 2.4
    distance = 4000
    cluster_mass = 10**5.
    deltaAKs = 0.05
    mass_sampling=5

    # Test filters
    filt_list = ['nirc2,J', 'nirc2,Kp']
    
    startTime = time.time()
    
    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere

    red_law = reddening.RedLawNishiyama09()
    
    iso = syn.IsochronePhot(logAge, AKs, distance,
                            evo_model=evo, atm_func=atm_func,
                            red_law=red_law, filters=filt_list,
                                mass_sampling=mass_sampling)

    print('Constructed isochrone: %d seconds' % (time.time() - startTime))

    imf_mass_limits = np.array([0.07, 0.5, 1, np.inf])
    imf_powers = np.array([-1.3, -2.3, -2.3])

    ##########
    # Start without multiplicity
    ##########
    my_imf1 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=None)
    print('Constructed IMF: %d seconds' % (time.time() - startTime))
    
    cluster1 = syn.ResolvedClusterDiffRedden(iso, my_imf1, cluster_mass, deltaAKs)
    clust1 = cluster1.star_systems
    print('Constructed cluster: %d seconds' % (time.time() - startTime))

    assert len(clust1) > 0
    
    plt.figure(3)
    plt.clf()
    plt.plot(clust1['m_nirc2_J'] - clust1['m_nirc2_Kp'], clust1['m_nirc2_J'], 'r.')
    plt.plot(iso.points['m_nirc2_J'] - iso.points['m_nirc2_Kp'], iso.points['m_nirc2_J'], 'c.')
    plt.gca().invert_yaxis()

    # *** Visual Inspections: ***
    #  - check that points (red) fall between isochrone points (blue)

    ##########
    # Test with multiplicity
    ##########
    multi = multiplicity.MultiplicityUnresolved()
    my_imf2 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=multi)
    print('Constructed IMF with multiples: %d seconds' % (time.time() - startTime))
    
    cluster2 = syn.ResolvedClusterDiffRedden(iso, my_imf2, cluster_mass, deltaAKs)
    clust2 = cluster2.star_systems
    print('Constructed cluster with multiples: %d seconds' % (time.time() - startTime))

    assert len(clust2) > 0
    assert len(cluster2.companions) > 0
    assert np.sum(clust2['N_companions']) == len(cluster2.companions)

    ##########
    # Plots 
    ##########
    # Plot an IR CMD and compare cluster members to isochrone.
    plt.figure(1)
    plt.clf()
    plt.plot(clust1['m_nirc2_J'] - clust1['m_nirc2_Kp'], clust1['m_nirc2_J'], 'r.')
    plt.plot(clust2['m_nirc2_J'] - clust2['m_nirc2_Kp'], clust2['m_nirc2_J'], 'b.')
    plt.plot(iso.points['m_nirc2_J'] - iso.points['m_nirc2_Kp'], iso.points['m_nirc2_J'], 'c-')
    plt.gca().invert_yaxis()
    plt.xlabel('J - Kp (mag)')
    plt.ylabel('J (mag')

    # Plot a mass-magnitude relationship.
    plt.figure(2)
    plt.clf()
    plt.semilogx(clust1['mass'], clust1['m_nirc2_J'], 'r.')
    plt.semilogx(clust2['mass'], clust2['m_nirc2_J'], 'r.')
    plt.gca().invert_yaxis()
    plt.xlabel('Mass (Msun)')
    plt.ylabel('J (mag)')

    return
    
def test_UnresolvedCluster():
    log_age = 6.7
    AKs = 0.0
    distance = 4000
    metallicity=0
    cluster_mass = 10**4.

    startTime = time.time()    
    multi = multiplicity.MultiplicityUnresolved()
    imf_in = imf.Kroupa_2001(multiplicity=multi)
    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere
    iso = syn.Isochrone(log_age, AKs, distance, metallicity=metallicity,
                            evo_model=evo, atm_func=atm_func, mass_sampling=10)
    print('Made Isochrone: %d seconds' % (time.time() - startTime))

    cluster = syn.UnresolvedCluster(iso, imf_in, cluster_mass)
    print('Constructed unresolved cluster: %d seconds' % (time.time() - startTime))

    # Plot an integrated spectrum of the whole cluster.
    wave = cluster.wave_trim
    flux = cluster.spec_trim
    plt.clf()
    plt.plot(wave, flux, 'k.')

    return

def test_ifmr_multiplicity():
    # Define cluster parameters
    logAge = 9.7
    AKs = 0.0
    distance = 1000
    cluster_mass = 1e6
    mass_sampling = 5

    # Test all filters
    filt_list = ['nirc2,Kp', 'nirc2,H', 'nirc2,J']

    startTime = time.time()
    
    evo = evolution.MISTv1()
    atm_func = atmospheres.get_merged_atmosphere
    ifmr_obj = ifmr.IFMR_Raithel18()

    red_law = reddening.RedLawNishiyama09()
    
    iso = syn.IsochronePhot(logAge, AKs, distance,
                            evo_model=evo, atm_func=atm_func,
                            red_law=red_law, filters=filt_list,
                            mass_sampling=mass_sampling)

    print('Constructed isochrone: %d seconds' % (time.time() - startTime))

    # Now to create the cluster.
    imf_mass_limits = np.array([0.07, 0.5, 1, np.inf])
    imf_powers = np.array([-1.3, -2.3, -2.3])

    ##########
    # Start without multiplicity and IFMR
    ##########
    my_imf1 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=None)
    print('Constructed IMF: %d seconds' % (time.time() - startTime)) 
    
    cluster1 = syn.ResolvedCluster(iso, my_imf1, cluster_mass, ifmr=ifmr_obj)
    clust1 = cluster1.star_systems
    print('Constructed cluster: %d seconds' % (time.time() - startTime))

   
    ##########
    # Test with multiplicity and IFMR
    ##########
    multi = multiplicity.MultiplicityUnresolved()
    my_imf2 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=multi)
    print('Constructed IMF with multiples: %d seconds' % (time.time() - startTime))
    
    cluster2 = syn.ResolvedCluster(iso, my_imf2, cluster_mass, ifmr=ifmr_obj)
    clust2 = cluster2.star_systems
    comps2 = cluster2.companions
    print('Constructed cluster with multiples: %d seconds' % (time.time() - startTime))

    ##########
    # Tests
    ##########

    # Check that we have black holes, neutron stars, and white dwarfs in both.
    assert len(np.where(clust1['phase'] == 101)) > 0   # WD
    assert len(np.where(clust2['phase'] == 101)) > 0
    assert len(np.where(clust1['phase'] == 102)) > 0   # NS
    assert len(np.where(clust2['phase'] == 102)) > 0
    assert len(np.where(clust1['phase'] == 103)) > 0   # BH
    assert len(np.where(clust2['phase'] == 103)) > 0

    # Now check that we have companions that are WDs, NSs, and BHs
    assert len(np.where(comps2['phase'] == 101)) > 0
    assert len(np.where(comps2['phase'] == 102)) > 0
    assert len(np.where(comps2['phase'] == 103)) > 0

    # Make sure no funky phase designations (due to interpolation effects)
    # slipped through
    idx = np.where( (clust1['phase'] > 5) & (clust1['phase'] < 101) & (clust1['phase'] != 9) )
    idx2 = np.where( (comps2['phase'] > 5) & (comps2['phase'] < 101) & (comps2['phase'] != 9) )
    assert len(idx[0]) == 0

    return

def test_metallicity():
    """
    Test isochrone generation at different metallicities
    """
    # Define isochrone parameters
    logAge = np.log10(5*10**6.) 
    AKs = 0.8 
    dist = 4000 
    evo_model = evolution.MISTv1()
    atm_func = atmospheres.get_phoenixv16_atmosphere
    red_law = reddening.RedLawHosek18b()
    filt_list = ['wfc3,ir,f127m', 'wfc3,ir,f139m', 'wfc3,ir,f153m']

    # Start with a solar metallicity isochrone    
    metallicity= 0.0

    # Make Isochrone object, with high mass_sampling to decrease compute time
    my_iso = syn.IsochronePhot(logAge, AKs, dist, metallicity=metallicity,
                            evo_model=evo_model, atm_func=atm_func,
                            red_law=red_law, filters=filt_list,
                            mass_sampling=10)

    # Test isochrone properties
    assert my_iso.points.meta['METAL_IN'] == 0.0
    assert os.path.exists('iso_6.70_0.80_04000_p00.fits')

    # Now for non-solar metallicity
    metallicity= -1.5

    # Make Isochrone object, with high mass_sampling to decrease compute time
    my_iso = syn.IsochronePhot(logAge, AKs, dist, metallicity=metallicity,
                            evo_model=evo_model, atm_func=atm_func,
                            red_law=red_law, filters=filt_list,
                            mass_sampling=10)

    # MIST model sub-directory names changed in SPISEA v2.1.4 update;
    # changing what "metal_act" value was. version 1 of MIST grid
    # names mistakenly used 0.015 for Zsolar, which was fixed for
    # version 2 onward. Note that this only effects this variable,
    # and not the grid output itself.
    # So, here the expected value for metal_act changes depending
    # on the grid value
    grid_version = evolution.get_installed_grid_num('{0}/evolution/'.format(os.environ['SPISEA_MODELS']))

    if grid_version == 1:
        metal_act = np.log10(0.00047 / 0.0142) # For Mist isochrones
    else:
        metal_act = np.log10(0.00045 / 0.0142) # For Mist isochrones
        
    # Test isochrone properties
    assert my_iso.points.meta['METAL_IN'] == -1.5
    assert np.isclose(my_iso.points.meta['METAL_ACT'], metal_act)
    assert os.path.exists('iso_6.70_0.80_04000_m15.fits')
    
    return

def test_cluster_mass():
    # Define cluster parameters
    logAge = 6.7
    AKs = 2.4
    distance = 4000
    cluster_mass = 10**5.
    mass_sampling = 5

    # Test filters
    filt_list = ['nirc2,J', 'nirc2,Kp']

    startTime = time.time()
    
    # Define evolution/atmosphere models and extinction law
    evo = evolution.MISTv1() 
    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawHosek18b()
    
    iso = syn.IsochronePhot(logAge, AKs, distance,
                            evo_model=evo, atm_func=atm_func,
                            red_law=red_law, filters=filt_list,
                            mass_sampling=mass_sampling)

    print('Constructed isochrone: %d seconds' % (time.time() - startTime))

    # Now to create the cluster.
    imf_mass_limits = np.array([0.2, 0.5, 1, 120.0])
    imf_powers = np.array([-1.3, -2.3, -2.3])

    # IFMR
    my_ifmr = ifmr.IFMR_Raithel18()
    

    ##########
    # Start without multiplicity
    ##########
    my_imf1 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=None)
    print('Constructed IMF: %d seconds' % (time.time() - startTime))
    
    cluster1 = syn.ResolvedCluster(iso, my_imf1, cluster_mass, ifmr=my_ifmr)
    clust1 = cluster1.star_systems
    print('Constructed cluster: %d seconds' % (time.time() - startTime))

    # Check that the total mass is within tolerance of input mass
    cluster_mass_out = clust1['systemMass'].sum()
    assert np.abs(cluster_mass_out - cluster_mass) < 200.0   # within 200 Msun of desired mass.
    print('Cluster Mass: IN = ', cluster_mass, " OUT = ", cluster_mass_out)

    ##########
    # Test with multiplicity
    ##########
    multi = multiplicity.MultiplicityUnresolved()
    my_imf2 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=multi)
    print('Constructed IMF with multiples: %d seconds' % (time.time() - startTime))
    
    cluster2 = syn.ResolvedCluster(iso, my_imf2, cluster_mass, ifmr=my_ifmr)
    clust2 = cluster2.star_systems
    print('Constructed cluster with multiples: %d seconds' % (time.time() - startTime))

    # Check that the total mass is within tolerance of input mass
    cluster_mass_out = clust2['systemMass'].sum()
    assert np.abs(cluster_mass_out - cluster_mass) < 200.0   # within 200 Msun of desired mass.
    print('Cluster Mass: IN = ', cluster_mass, " OUT = ", cluster_mass_out)

    return

def test_compact_object_companions():
    
    # Define cluster parameters
    logAge = 6.7
    AKs = 2.4
    distance = 4000
    cluster_mass = 10**4.
    mass_sampling=5

    # Test filters
    filt_list = ['nirc2,J', 'nirc2,Kp']

    startTime = time.time()
    
    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere

    red_law = reddening.RedLawNishiyama09()
    
    iso = syn.IsochronePhot(logAge, AKs, distance,
                            evo_model=evo, atm_func=atm_func,
                            red_law=red_law, filters=filt_list,
                            mass_sampling=mass_sampling)

    print('Constructed isochrone: %d seconds' % (time.time() - startTime))
    
    clust_multiplicity = multiplicity.MultiplicityResolvedDK()

    massLimits = np.array([0.2, 0.5, 1, 120]) # mass segments
    powers = np.array([-1.3, -2.3, -2.3]) # power-law exponents
    clust_imf_Mult = imf.IMF_broken_powerlaw(massLimits, powers, multiplicity=clust_multiplicity)
    #clust_imf_Mult = imf.Kroupa_2001(multiplicity=clust_multiplicity)
    clust_Mult = syn.ResolvedCluster(iso, clust_imf_Mult, cluster_mass, ifmr=ifmr.IFMR_Raithel18())

    # Makes sure compact object companions not including MIST WDs
    #(i.e. those with no luminosity) are being given phases
    nan_lum_companions = clust_Mult.companions[np.isnan(clust_Mult.companions['L'])]

    assert (len(nan_lum_companions) == 0) | (all(np.isnan(nan_lum_companions['phase'])) == False)

#=================================#
# Additional timing functions
#=================================#

def time_test_cluster():
    logAge = 6.7
    AKs = 2.7
    distance = 4000
    cluster_mass = 10**4

    startTime = time.time()
    
    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawNishiyama09()
    filt_list = ['nirc2,J', 'nirc2,Kp']
    
    iso = syn.IsochronePhot(logAge, AKs, distance,
                            evo_model=evo, atm_func=atm_func,
                            red_law=red_law, filters=filt_list)
    print('Constructed isochrone: %d seconds' % (time.time() - startTime))

    imf_limits = np.array([0.07, 0.5, 150])
    imf_powers = np.array([-1.3, -2.35])
    multi = multiplicity.MultiplicityUnresolved()
    my_imf = imf.IMF_broken_powerlaw(imf_limits, imf_powers, multiplicity=multi)
    print('Constructed IMF with multiples: %d seconds' % (time.time() - startTime))
    
    cluster = syn.ResolvedCluster(iso, my_imf, cluster_mass)
    print('Constructed cluster: %d seconds' % (time.time() - startTime))

    return
    
def model_young_cluster_object(resolved=False):
    log_age = 6.5
    AKs = 0.1
    distance = 8000.0
    cluster_mass = 10000.
    
    multi = multiplicity.MultiplicityUnresolved()
    imf_in = imf.Kroupa_2001(multiplicity=multi)
    evo = evolution.MergedPisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere
    iso = syn.Isochrone(log_age, AKs, distance, evo, mass_sampling=10)

    if resolved:
        cluster = syn.ResolvedCluster(iso, imf_in, cluster_mass)
    else:
        cluster = syn.UnresolvedCluster(iso, imf_in, cluster_mass, wave_range=[19000,24000])

    # Plot the spectrum of the most massive star
    idx = cluster.mass_all.argmax()
    print('Most massive star is {0:f} M_sun.'.format(cluster.mass_all[idx]))
    #bigstar = cluster.spec_list_trim[idx]
    plt.figure(1)
    plt.clf()
    plt.plot(cluster.spec_list_trim[idx]._wavetable, cluster.spec_list_trim[idx]._fluxtable, 'k.')

    # Plot an integrated spectrum of the whole cluster.
    wave, flux = cluster.spec_list_trim[idx]._wavetable, cluster.spec_trim
    plt.figure(2)
    plt.clf()
    plt.plot(wave, flux, 'k.')

    return
    
def time_test_mass_match():
    log_age = 6.7
    AKs = 2.7
    distance = 4000
    cluster_mass = 5e3
    
    imf_in = imf.Kroupa_2001(multiplicity=None)

    start_time = time.time()
    iso = syn.IsochronePhot(log_age, AKs, distance)
    iso_masses = iso.points['mass']
    print('Generated iso masses in {0:.0f} s'.format(time.time() - start_time))

    start_time = time.time()
    star_masses, isMulti, compMass, sysMass = imf_in.generate_cluster(cluster_mass)
    print('Generated cluster masses in {0:.0f} s'.format(time.time() - start_time))
    
    def match_model_masses1(isoMasses, starMasses):
        indices = np.empty(len(starMasses), dtype=int)
        
        for ii in range(len(starMasses)):
            theMass = starMasses[ii]
            
            dm = np.abs(isoMasses - theMass)
            mdx = dm.argmin()

            # Model mass has to be within 10% of the desired mass
            if (dm[mdx] / theMass) > 0.1:
                indices[ii] = -1
            else:
                indices[ii] = mdx

        return indices
            

    def match_model_masses2(isoMasses, starMasses):
        isoMasses_tmp = isoMasses.reshape((len(isoMasses), 1))
        kdt = KDTree(isoMasses_tmp)
 
        starMasses_tmp = starMasses.reshape((len(starMasses), 1))
        q_results = kdt.query(starMasses_tmp, k=1)
        indices = q_results[1]

        dm_frac = np.abs(starMasses - isoMasses[indices]) / starMasses

        idx = np.where(dm_frac > 0.1)[0]
        indices[idx] = -1
        
        return indices

    print('Test #1 START')
    start_time = time.time()
    idx1 = match_model_masses1(iso_masses, star_masses)
    stop_time = time.time()
    print('Test #1 STOPPED after {0:.0f} seconds'.format(stop_time - start_time))

    print('Test #2 START')
    start_time = time.time()
    idx2 = match_model_masses2(iso_masses, star_masses)
    stop_time = time.time()
    print('Test #2 STOPPED after {0:.0f} seconds'.format(stop_time - start_time))

    return

def generate_Spera15_IFMR():
    """
    Make a set of objects based on the Spera15 IFMR for the purposes of testing
    Expect 28 total objects
    8 invalids, 8 WDs, 2 NSs, and 10 BHs
    """
    Spera = ifmr.IFMR_Spera15()

    def FeH_from_Z(Z):
        return np.log10(Z/0.019)

    metal = np.array([2.0e-4, 1.0e-3, 2.0e-3, 2.0e-2]) #ensure that all Spera metallicity regimes are represented

    FeH = FeH_from_Z(metal) #generate death mass takes metallicty as [Fe/H]
    #want to get a good range of masses for Spera, should expect 8 invalids, 8 WDs, 3 NSs, and 9 BHs 
    ZAMS = np.array([-0.2*np.ones(len(FeH)), 0.2*np.ones(len(FeH)), 4.0*np.ones(len(FeH)), 9.2*np.ones(len(FeH)),
                    15.0*np.ones(len(FeH)), 30.0*np.ones(len(FeH)), 150.0*np.ones(len(FeH))])


    output_array = np.concatenate((Spera.generate_death_mass(ZAMS[0], FeH), Spera.generate_death_mass(ZAMS[1], FeH),
                                  Spera.generate_death_mass(ZAMS[2], FeH), Spera.generate_death_mass(ZAMS[3], FeH),
                                  Spera.generate_death_mass(ZAMS[4], FeH), Spera.generate_death_mass(ZAMS[5], FeH),
                                  Spera.generate_death_mass(ZAMS[6], FeH)), axis=1)

    #count up number of objects formed and ensure it matahces the number of stars input
    bad_idx = np.where(output_array[1] == -1)
    WD_idx = np.where(output_array[1] == 101)
    NS_idx = np.where(output_array[1] == 102)
    BH_idx = np.where(output_array[1] == 103)

    rem_mass = output_array[0]

    return bad_idx[0], WD_idx[0], NS_idx[0], BH_idx[0], rem_mass

def test_Spera15_IFMR1():
    """
    Check to make sure the total number of objects input matches the number of objects output (28)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    total = len(WD_idx) + len(NS_idx) + len(BH_idx) + len(bad_idx)
    assert total == 28 , "The # of input objects does not match the number of output objects for the Spera15 IFMR"

    return

def test_Spera15_IFMR2():
    """
    Check that all negative remnant masses have type code -1
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    #check to make sure no unhandled negative remnant masses (ie without type code -1 assigned)
    neg_mass_idx = np.where(rem_mass < 0)
    assert set(bad_idx) == set(neg_mass_idx[0]), "There are unhandled negative remnant masses for the Spera15 IFMR"

    return

def test_Spera15_IFMR_3():
    """
    Check that there are no left over zeroes in the remnant mass array
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    assert np.all(rem_mass != 0) , "There are left over zeros in the remnant mass array for the Spera15 IFMR"

    return

def test_Spera15_IFMR_4():
    """
    Check that the correct number of invalid objects were generated (8)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    assert len(bad_idx) == 8 , "There are not the right number of invalid objects for the Spera15 IFMR"

    return

def test_Spera15_IFMR_5():
    """
    Check that the correct number of WDs were generated (8)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    assert len(WD_idx) == 8 , "There are not the right number of WDs for the Spera15 IFMR"

    return

def test_Spera15_IFMR_6():
    """
    Check that the correct number of NSs were generated (8)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    assert len(NS_idx) == 2 , "There are not the right number of NSs for the Spera15 IFMR"

    return

def test_Spera15_IFMR_7():
    """
    Check that the correct number of invalid objects were generated (8)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    assert len(BH_idx) == 10 , "There are not the right number of BHs for the Spera15 IFMR"

    return
    
def generate_Raithel18_IFMR():
    """
    Make a set of objects using the Raithel18 IFMR for the purposes of testing
    Will make a total of 16 objects
    Will have 3 invalid objects and 2 WDs
    """

    Raithel = ifmr.IFMR_Raithel18()
    ZAMS = np.array([-0.2, 0.2, 1.0, 7.0, 10.0, 14.0, 16.0, 18.0, 18.6, 22.0, 26.0, 28.0, 50.0, 61.0, 119.0, 121.0]) 
    #3 invalid indices, 2 WDs, cannot make statements about #of BHs and NSs because the Raithel IFMR has some randomness

    output_array = Raithel.generate_death_mass(ZAMS)

    #count up number of objects formed and ensure it matahces the number of stars input
    bad_idx = np.where(output_array[1] == -1)
    WD_idx = np.where(output_array[1] == 101)
    NS_idx = np.where(output_array[1] == 102)
    BH_idx = np.where(output_array[1] == 103)

    rem_mass = output_array[0]

    return bad_idx[0], WD_idx[0], NS_idx[0], BH_idx[0], rem_mass

def test_Raithel18_IFMR_1():
    """
    Check to make sure that the number of input objects matches the number of output objects
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Raithel18_IFMR()

    total = len(WD_idx) + len(NS_idx) + len(BH_idx) + len(bad_idx)
    assert total == 16 , "The # of input objects does not match the number of output objects for the Raithel18 IFMR"

    return

def test_Raithel18_IFMR_2():
    """
    Check that all negative remnant masses have the type code -1
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Raithel18_IFMR()

    #check to make sure no unhandled negative remnant masses (ie without type code -1 assigned)
    neg_mass_idx = np.where(rem_mass < 0)
    assert set(bad_idx) == set(neg_mass_idx[0]) , "There are unhandled negative remnant masses for the Raithel18 IFMR"

    return

def test_Raithel18_IFMR_3():
    """
    Check that there are no left over zeros in the remnant mass array
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Raithel18_IFMR()

    #check to make sure there are no left over zeroes in the remnant mass array
    assert np.all(rem_mass != 0), "There are left over zeros in the remnant mass array for the Raithel18 IFMR"

    return

def test_Raithel18_IFMR_4():
    """
    Check that the right number of invalid objects are returned (3)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Raithel18_IFMR()

    assert len(bad_idx) == 3 , "There are not the right number of invalid objects for the Raithel18 IFMR"

    return

def test_Raithel18_IFMR_5():
    """
    Check that the right number of WDs are returned (2)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Raithel18_IFMR()

    assert len(WD_idx) == 2 , "There are not the right number of WDs for the Raithel18 IFMR"

    return
