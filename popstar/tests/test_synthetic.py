import time
import pylab as plt
import numpy as np
from popstar import synthetic
import pysynphot
from pysynphot import observation as obs
from pysynphot import spectrum
import pdb
from scipy.spatial import cKDTree as KDTree

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
    cluster_mass = 100000

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

    plt.figure(3)
    plt.clf()
    plt.plot(clust1['magJ'] - clust1['magKp'], clust1['magJ'], 'r.')
    plt.plot(iso.points['magJ'] - iso.points['magKp'], iso.points['magJ'], 'c.')
    plt.gca().invert_yaxis()

    # *** Visual Inspections: ***
    #  - check that points (red) fall between isochrone points (blue)

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
    from popstar import synthetic as syn
    from popstar import atmospheres as atm
    from popstar import evolution
    from popstar import reddening
    from popstar.imf import imf
    from popstar.imf import multiplicity

    logAge = 6.7
    AKs = 2.4
    distance = 4000
    cluster_mass = 100000
    deltaAKs = 0.05
    
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
    
    cluster1 = syn.ResolvedClusterDiffRedden(iso, my_imf1, cluster_mass, deltaAKs)
    clust1 = cluster1.star_systems
    print 'Constructed cluster: %d seconds' % (time.time() - startTime)

    plt.figure(3)
    plt.clf()
    plt.plot(clust1['magJ'] - clust1['magKp'], clust1['magJ'], 'r.')
    plt.plot(iso.points['magJ'] - iso.points['magKp'], iso.points['magJ'], 'c.')
    plt.gca().invert_yaxis()

    # *** Visual Inspections: ***
    #  - check that points (red) fall between isochrone points (blue)

    ##########
    # Test with multiplicity
    ##########
    multi = multiplicity.MultiplicityUnresolved()
    my_imf2 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=multi)
    print 'Constructed IMF with multiples: %d seconds' % (time.time() - startTime)
    
    cluster2 = syn.ResolvedClusterDiffRedden(iso, my_imf2, cluster_mass, deltaAKs)
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
    evo = evolution.MergedBaraffePisaEkstromParsec()
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

def time_test_cluster():
    from popstar import synthetic as syn
    from popstar import atmospheres as atm
    from popstar import evolution
    from popstar import reddening
    from popstar.imf import imf
    from popstar.imf import multiplicity

    logAge = 6.7
    AKs = 2.7
    distance = 4000
    cluster_mass = 5e3

    startTime = time.time()
    
    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atm.get_merged_atmosphere
    red_law = reddening.RedLawNishiyama09()
    
    iso = syn.IsochronePhot(logAge, AKs, distance,
                            evo_model=evo, atm_func=atm_func,
                            red_law=red_law)
    print 'Constructed isochrone: %d seconds' % (time.time() - startTime)

    imf_limits = np.array([0.07, 0.5, 150])
    imf_powers = np.array([-1.3, -2.35])
    multi = multiplicity.MultiplicityUnresolved()
    my_imf = imf.IMF_broken_powerlaw(imf_limits, imf_powers, multiplicity=multi)
    print 'Constructed IMF with multiples: %d seconds' % (time.time() - startTime)
    
    cluster = syn.ResolvedCluster(iso, my_imf, cluster_mass)
    print 'Constructed cluster: %d seconds' % (time.time() - startTime)

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
    
def time_test_mass_match():
    from popstar import synthetic as syn
    from popstar import atmospheres as atm
    from popstar import evolution
    from popstar.imf import imf
    from popstar.imf import multiplicity

    log_age = 6.7
    AKs = 2.7
    distance = 4000
    cluster_mass = 5e3
    
    imf_in = imf.Kroupa_2001(multiplicity=None)

    start_time = time.time()
    iso = syn.IsochronePhot(log_age, AKs, distance)
    iso_masses = iso.points['mass']
    print 'Generated iso masses in {0:.0f} s'.format(time.time() - start_time)

    start_time = time.time()
    star_masses, isMulti, compMass, sysMass = imf_in.generate_cluster(cluster_mass)
    print 'Generated cluster masses in {0:.0f} s'.format(time.time() - start_time)
    
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

    print 'Test #1 START'
    start_time = time.time()
    idx1 = match_model_masses1(iso_masses, star_masses)
    stop_time = time.time()
    print 'Test #1 STOPPED after {0:.0f} seconds'.format(stop_time - start_time)

    print 'Test #2 START'
    start_time = time.time()
    idx2 = match_model_masses2(iso_masses, star_masses)
    stop_time = time.time()
    print 'Test #2 STOPPED after {0:.0f} seconds'.format(stop_time - start_time)

    return

def test_binning_methods():
    """
    Compare the pysynphot binflux routine on filter integrations at different resolutions.
    Shows bug in routine: integrated flux depends on the filter resolution! 

    Fix: manually integrate the binned filter function. Shows that this method performs
    much better, getting nearly the same output flux for different filter resolutions,
    as we would expect
    """
    # We'll test an integration of the vega spectrum through the WFC3-IR F127M filter
    vega = synthetic.Vega()
    filt = pysynphot.ObsBandpass('wfc3,ir,f127m')

    # Convert to ArraySpectralElement for resampling.
    filt = spectrum.ArraySpectralElement(filt.wave, filt.throughput,
                                             waveunits=filt.waveunits)
    
    # Two rebinning schemes: one coarse and the other fine
    idx = np.where(filt.throughput > 0.001)[0]
    new_wave = np.linspace(filt.wave[idx[0]], filt.wave[idx[-1]], 1500, dtype=float)
    filt_fine = filt.resample(new_wave)

    wave_bin = vega.wave
    filt_bin = synthetic.rebin_spec(filt.wave, filt.throughput, wave_bin)        
    filt_coarse = pysynphot.ArrayBandpass(wave_bin, filt_bin)

    # Do the filter integration in 2 methods: one with pysynphot binflux,
    # the other with manual integration
    vega_obs_fine = obs.Observation(vega, filt_fine, binset=filt_fine.wave, force='taper')
    vega_obs_coarse = obs.Observation(vega, filt_coarse, binset=filt_coarse.wave, force='taper')

    fine_binflux = vega_obs_fine.binflux.sum()
    coarse_binflux = vega_obs_coarse.binflux.sum()

    diff_f = np.diff(vega_obs_fine.binwave)
    diff_f = np.append(diff_f, diff_f[-1])
    fine_manual = np.sum(vega_obs_fine.binflux * diff_f)
    
    diff_c = np.diff(vega_obs_coarse.binwave)
    diff_c = np.append(diff_c, diff_c[-1])
    coarse_manual = np.sum(vega_obs_coarse.binflux * diff_c)

    print('**************************************')
    print('Integrated flux with binflux:')
    print('fine binning: {0}'.format(fine_binflux))
    print('coarse binning: {0}'.format(coarse_binflux))
    print('And with manual integration:')
    print('fine binning: {0}'.format(fine_manual))
    print('coarse binning: {0}'.format(coarse_manual))
    print('**************************************')

    pdb.set_trace()
    return
