import time
import pylab as plt
import numpy as np
import pdb

def test_isochrone(plot=False):
    from popstar import synthetic as syn

    logAge = 6.7
    AKs = 2.7
    distance = 4000

    startTime = time.time()
    iso = syn.Isochrone(logAge, AKs, distance)
    print 'Test completed in: %d seconds' % (time.time() - startTime)
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

def test_IsochronePhot(plot=False):
    from popstar import synthetic as syn

    logAge = 6.7
    AKs = 2.7
    distance = 4000

    startTime = time.time()
    iso = syn.IsochronePhot(logAge, AKs, distance)
    endTime = time.time()
    print 'Test completed in: %d seconds' % (endTime - startTime)
    # Typically takes 120 seconds if file is regenerated.
    # Limited by pysynphot.Icat call in atmospheres.py

    assert iso.points.meta['LOGAGE'] == logAge
    assert iso.points.meta['AKS'] == AKs
    assert iso.points.meta['DISTANCE'] == distance
    assert len(iso.points) > 100

    assert 'magJ' in iso.points.colnames

    if plot:
        plt.figure(1) 
        iso.plot_CMD('mag814w', 'mag160w')
        
        plt.figure(2)
        iso.plot_mass_magnitude('mag160w')

    return

def test_ResolvedCluster():
    from popstar import synthetic as syn
    from popstar import atmospheres as atm
    from popstar import evolution
    from popstar import reddening
    from popstar.imf import imf
    from popstar.imf import multiplicity

    logAge = 6.7
    AKs = 2.4
    distance = 4000
    cluster_mass = 5000

    startTime = time.time()
    
    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atm.get_merged_atmosphere

    red_law = reddening.RedLawNishiyama09()
    
    iso = syn.IsochronePhot(logAge, AKs, distance,
                            evo_model=evo, atm_func=atm_func,
                            red_law=red_law)

    print 'Constructed isochrone: %d seconds' % (time.time() - startTime)

    imf_mass_limits = np.array([0.07, 0.5, 1, np.inf])
    imf_powers = np.array([-1.3, -2.3, -2.3])

    ##########
    # Start without multiplicity
    ##########
    my_imf1 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=None)
    print 'Constructed IMF: %d seconds' % (time.time() - startTime)
    
    cluster1 = syn.ResolvedCluster(iso, my_imf1, cluster_mass)
    clust1 = cluster1.star_systems
    print 'Constructed cluster: %d seconds' % (time.time() - startTime)

    plt.clf()
    plt.plot(clust1['magJ'] - clust1['magKp'], clust1['magJ'], 'r.')
    plt.plot(iso.points['magJ'] - iso.points['magKp'], iso.points['magJ'], 'c.')
    plt.gca().invert_yaxis()

    ##########
    # Test with multiplicity
    ##########
    multi = multiplicity.MultiplicityUnresolved()
    my_imf2 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=multi)
    print 'Constructed IMF with multiples: %d seconds' % (time.time() - startTime)
    
    cluster2 = syn.ResolvedCluster(iso, my_imf2, cluster_mass)
    clust2 = cluster2.star_systems
    print 'Constructed cluster with multiples: %d seconds' % (time.time() - startTime)

    ##########
    # Plots 
    ##########
    # Plot an IR CMD and compare cluster members to isochrone.
    plt.figure(1)
    plt.clf()
    plt.plot(clust1['magJ'] - clust1['magKp'], clust1['magJ'], 'r.')
    plt.plot(clust2['magJ'] - clust2['magKp'], clust2['magJ'], 'b.')
    plt.plot(iso.points['magJ'] - iso.points['magKp'], iso.points['magJ'], 'c-')
    plt.gca().invert_yaxis()
    plt.xlabel('J - Kp (mag)')
    plt.ylabel('J (mag')

    # Plot a mass-magnitude relationship.
    plt.figure(2)
    plt.clf()
    plt.semilogx(clust1['mass'], clust1['magJ'], 'r.')
    plt.semilogx(clust2['mass'], clust2['magJ'], 'r.')
    plt.gca().invert_yaxis()
    plt.xlabel('Mass (Msun)')
    plt.ylabel('J (mag)')
    
    return clust2

    # # Plot the spectrum of the most massive star
    # idx = cluster.mass.argmax()
    # plt.clf()
    # plt.plot(cluster.stars[idx].wave, cluster.stars[idx].flux, 'k.')

    # # Plot an integrated spectrum of the whole cluster.
    # wave, flux = cluster.get_integrated_spectrum()
    # plt.clf()
    # plt.plot(wave, flux, 'k.')

    return
    
def test_UnresolvedCluster():
    from popstar import synthetic as syn
    from popstar import atmospheres as atm
    from popstar import evolution
    from popstar.imf import imf
    from popstar.imf import multiplicity
    
    log_age = 6.7
    AKs = 0.0
    distance = 4000
    cluster_mass = 30.

    startTime = time.time()    
    multi = multiplicity.MultiplicityUnresolved()
    imf_in = imf.Kroupa_2001(multiplicity=multi)
    evo = evolution.MergedPisaEkstromParsec()
    iso = syn.Isochrone(log_age, AKs, distance, evo, mass_sampling=10)
    print 'Made cluster: %d seconds' % (time.time() - startTime)

    cluster = syn.UnresolvedCluster(iso, imf_in, cluster_mass)
    print 'Constructed unresolved cluster: %d seconds' % (time.time() - startTime)

    # Plot an integrated spectrum of the whole cluster.
    wave = cluster.spec_trim.wave
    flux = cluster.spec_trim.flux
    plt.clf()
    plt.plot(wave, flux, 'k.')
    pdb.set_trace()
    return

def model_young_cluster_object(resolved=False):
    from popstar import synthetic as syn
    from popstar import atmospheres as atm
    from popstar import evolution
    from popstar.imf import imf
    from popstar.imf import multiplicity

    log_age = 6.5
    AKs = 0.1
    distance = 8000.0
    cluster_mass = 10000.
    
    multi = multiplicity.MultiplicityUnresolved()
    imf_in = imf.Kroupa_2001(multiplicity=multi)
    evo = evolution.MergedPisaEkstromParsec()
    atm_func = atm.get_merged_atmosphere
    iso = syn.Isochrone(log_age, AKs, distance, evo, mass_sampling=10)

    if resolved:
        cluster = syn.ResolvedCluster(iso, imf_in, cluster_mass)
    else:
        cluster = syn.UnresolvedCluster(iso, imf_in, cluster_mass, wave_range=[19000,24000])

    # Plot the spectrum of the most massive star
    idx = cluster.mass_all.argmax()
    print 'Most massive star is {0:f} M_sun.'.format(cluster.mass_all[idx])
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

def time_stuff():
    time.time()

    return
    
