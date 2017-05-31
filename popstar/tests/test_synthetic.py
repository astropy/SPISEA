import time
import pylab as plt
import numpy as np
from popstar import synthetic, reddening
import pysynphot
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


def test_phot_consistency(filt='all'):
    """
    Test photometric consistency of generated isochrone (IsochronePhot)
    against pre-generated isochrone with native filter sampling. Requires
    consistency to within 0.005 mag.

    Base isochrone is at 5 Myr, AKs = 0, 1000 pc, mass_sampling=10 

    Paramters:
    ----------
    filt: 'all', 'hst', 'vista', 'decam', 'ps1', 'jwst'
        Specify what filter set you want to test

    """
    from astropy.table import Table
    import os
    
    # Load pre-generated isochrone, located in popstar tests directory
    direct = os.path.dirname(__file__)
    orig = Table.read(direct+'/iso_6.70_0.00_01000.fits', format='fits')
    #orig = Table.read(direct+'/iso_6.40_2.40_08000.fits', format='fits')

    # Generate new isochrone with popstar code
    if filt == 'all':
        filt_list = {'hst_F127M': 'wfc3,ir,f127m', 'hst_F139M': 'wfc3,ir,f139m', 'hst_F153M': 'wfc3,ir,f153m',
                         'hst_F814W': 'acs,wfc1,f814w', 'hst_F125W': 'wfc3,ir,f125w', 'hst_F160W': 'wfc3,ir,f160w',
                         'decam_y': 'decam,y', 'decam_i': 'decam,i', 'decam_z': 'decam,z',
                         'decam_u':'decam,u', 'decam_g':'decam,g', 'decam_r':'decam,r',
                         'vista_Y':'vista,Y', 'vista_Z':'vista,Z', 'vista_J': 'vista,J',
                         'vista_H': 'vista,H', 'vista_Ks': 'vista,Ks',
                         'ps1_z':'ps1,z', 'ps1_g':'ps1,g', 'ps1_r': 'ps1,r',
                         'ps1_i': 'ps1,i', 'ps1_y':'ps1,y',
                         'jwst_F090W': 'jwst,F090W', 'jwst_F164N': 'jwst,F164N', 'jwst_F212N': 'jwst,F212N',
                         'jwst_F323N':'jwst,F323N', 'jwst_F466N': 'jwst,F466N',
                         'nirc2_J': 'nirc2,J', 'nirc2_H': 'nirc2,H', 'nirc2_Kp': 'nirc2,Kp', 'nirc2_K': 'nirc2,K',
                         'nirc2_Lp': 'nirc2,Lp', 'nirc2_Hcont': 'nirc2,Hcont',
                         'nirc2_FeII': 'nirc2,FeII', 'nirc2_Brgamma': 'nirc2,Brgamma'}

    elif filt == 'decam':
        filt_list = {'decam_y': 'decam,y', 'decam_i': 'decam,i', 'decam_z': 'decam,z',
                         'decam_u':'decam,u', 'decam_g':'decam,g', 'decam_r':'decam,r'}
    elif filt == 'vista':
        filt_list = {'vista_Y':'vista,Y', 'vista_Z':'vista,Z', 'vista_J': 'vista,J',
                         'vista_H': 'vista,H', 'vista_Ks': 'vista,Ks'}
    elif filt == 'ps1':
        filt_list = {'ps1_z':'ps1,z', 'ps1_g':'ps1,g', 'ps1_r': 'ps1,r', 'ps1_i': 'ps1,i',
                         'ps1_y':'ps1,y'}
    elif filt == 'jwst':
        filt_list = {'jwst_F090W': 'jwst,F090W', 'jwst_F164N': 'jwst,F164N', 'jwst_F212N': 'jwst,F212N',
                         'jwst_F323N':'jwst,F323N', 'jwst_F466N': 'jwst,F466N'}
    elif filt == 'hst':
        filt_list = {'hst_F127M': 'wfc3,ir,f127m', 'hst_F139M': 'wfc3,ir,f139m', 'hst_F153M': 'wfc3,ir,f153m',
                         'hst_F814W': 'acs,wfc1,f814w', 'hst_F125W': 'wfc3,ir,f125w', 'hst_F160W': 'wfc3,ir,f160w'}
    elif filt == 'nirc2':
        filt_list = {'nirc2_J': 'nirc2,J', 'nirc2_H': 'nirc2,H', 'nirc2_Kp': 'nirc2,Kp', 'nirc2_K': 'nirc2,K',
                         'nirc2_Lp': 'nirc2,Lp','nirc2_Ms': 'nirc2,Ms', 'nirc2_Hcont': 'nirc2,Hcont',
                         'nirc2_FeII': 'nirc2,FeII', 'nirc2_Brgamma': 'nirc2,Brgamma'}

    print 'Making isochrone'
    iso = synthetic.IsochronePhot(6.7, 0, 1000, mass_sampling=10, filters=filt_list, rebin=True)
    #redlaw = reddening.RedLawNishiyama09()
    #iso = synthetic.IsochronePhot(6.4, 2.4, 8000, mass_sampling=10, red_law=redlaw, filters=filt_list, rebin=True)
    iso = iso.points

    # First assert that the stellar masses are the same
    diff = abs(orig['mass'] - iso['mass'])
    assert np.sum(diff) == 0

    # Identify the photometry columns
    cols = iso.keys()
    idx = []
    for ii in range(len(cols)):
        if cols[ii].startswith('mag'):
            idx.append(ii)
    mag_cols = np.array(cols)[idx]

    # Test the consistency of each column with the original isochrone
    for ii in mag_cols:
        orig_mag = orig[ii]
        new_mag = iso[ii]

        np.testing.assert_allclose(orig_mag, new_mag, rtol=0.01, err_msg="{0} failed".format(ii))

        # Also report median abs difference
        diff = abs(orig_mag - new_mag)
        print('{0} median abs diff: {1}'.format(ii, np.median(diff)))


    print('Phot consistency test successful for {0}'.format(filt))
    return

