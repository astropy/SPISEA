import numpy as np
import time
    
def test_create_MultiplicityUnresolved():
    """
    Tests creating and accessing a MultiplicityUnresolved object.
    """
    from .. import multiplicity

    # All default parameters -- check their values
    mu1 = multiplicity.MultiplicityUnresolved()
    assert mu1.MF_amp == 0.44
    assert mu1.MF_pow == 0.51
    assert mu1.CSF_amp == 0.50
    assert mu1.CSF_pow == 0.45
    assert mu1.CSF_max == 3
    assert mu1.q_pow == -0.4
    assert mu1.q_min == 0.01

    # Test setting different parameters
    mu2 = multiplicity.MultiplicityUnresolved(MF_amp=0.4, 
                                             MF_power=0.4,
                                             CSF_amp=0.4, 
                                             CSF_power=0.4, 
                                             CSF_max=4,
                                             q_power=0.4, 
                                             q_min=0.04)
    assert mu2.MF_amp == 0.4 
    assert mu2.MF_pow == 0.4
    assert mu2.CSF_amp == 0.4 
    assert mu2.CSF_pow == 0.4 
    assert mu2.CSF_max == 4
    assert mu2.q_pow == 0.4 
    assert mu2.q_min == 0.04


def test_multiplicity_fraction():
    """
    Test creating a MultiplicityUnresolved object and getting
    the multiplicity fraction out.
    """    
    from spisea.imf import multiplicity
    
    # First set of multiplicity parameters
    mu1 = multiplicity.MultiplicityUnresolved()

    mf1_1 = mu1.multiplicity_fraction(1.0)
    np.testing.assert_almost_equal(mf1_1, 0.44, decimal=2)

    mf1_2 = mu1.multiplicity_fraction(10.0)
    np.testing.assert_almost_equal(mf1_2, 1.0, decimal=2)

    mf1_3 = mu1.multiplicity_fraction(0.1)
    np.testing.assert_almost_equal(mf1_3, 0.136, decimal=2)

    # Second set of multiplicity parameters
    mu2 = multiplicity.MultiplicityUnresolved(MF_amp=0.4, MF_power=0.4,
                                             CSF_amp=0.4, CSF_power=0.4, CSF_max=4,
                                             q_power=0.4, q_min=0.04)

    mf2_1 = mu1.multiplicity_fraction(1.0)
    np.testing.assert_almost_equal(mf2_1, 0.44, decimal=2)

    mf2_2 = mu1.multiplicity_fraction(10.0)
    np.testing.assert_almost_equal(mf2_2, 1.0, decimal=2)

    mf2_3 = mu1.multiplicity_fraction(0.1)
    np.testing.assert_almost_equal(mf2_3, 0.136, decimal=2)

    # Test brown dwarf mass fractions
    mf_bd1 = mu1.multiplicity_fraction(0.07)  # near upper BD limit
    mf_bd2 = mu1.multiplicity_fraction(0.04)  # mid BD
    mf_bd3 = mu1.multiplicity_fraction(0.01)  # lower BD limit
    assert np.isclose(mf_bd1, 0.16, atol=0.01)
    assert np.isclose(mf_bd2, 0.08, atol=0.01)
    assert np.isclose(mf_bd3, 0.0, atol=1e-6)


def test_multiplicity_fraction_array():
    """
    Test multiplicity_fraction() on the MultiplicityUnresolved object
    where the inputs and outputs are in array form.
    """
    from spisea.imf import multiplicity
    
    # First set of multiplicity parameters
    mu1 = multiplicity.MultiplicityUnresolved()

    mass_array = np.array([1.0, 10.0, 0.1, 0.07, 0.04, 0.01])
    mf_array = mu1.multiplicity_fraction(mass_array)

    # Stellar regime checks
    np.testing.assert_almost_equal(mf_array[0], 0.44, decimal=2)
    np.testing.assert_almost_equal(mf_array[1], 1.0, decimal=2)
    np.testing.assert_almost_equal(mf_array[2], 0.136, decimal=2)

    # BD regime checks
    # interpolation between values implies lower masses --> lower mf
    assert mf_array[3] < mf_array[2]   
    assert mf_array[4] <= mf_array[3]  
    assert mf_array[5] <= mf_array[4]  

    # Ensure mf stars within reasonable bound (upper limit is 0.2)
    assert np.all(mf_array[3:] >= 0.0)
    assert np.all(mf_array[3:] <= 0.2)
    
    
def test_companion_star_fraction():
    """
    Test the companion_star fraction on the MultiplicityUnresolved object.
    """
    from spisea.imf import multiplicity

    # First set of multiplicity parameters
    mu1 = multiplicity.MultiplicityUnresolved()

    csf1_1 = mu1.companion_star_fraction(1.0)
    np.testing.assert_almost_equal(csf1_1, 0.5, decimal=2)

    csf1_2 = mu1.companion_star_fraction(70.0)
    np.testing.assert_almost_equal(csf1_2, 3.0, decimal=2)

    csf1_3 = mu1.companion_star_fraction(0.1)
    np.testing.assert_almost_equal(csf1_3, 0.177, decimal=2)

    # Second set of multiplicity parameters
    mu2 = multiplicity.MultiplicityUnresolved(MF_amp=0.4, MF_power=0.4,
                                              CSF_amp=0.4, CSF_power=0.4, CSF_max=2,
                                              q_power=0.4, q_min=0.04)

    # csf2_1 = mu1.companion_star_fraction(1.0)
    # np.testing.assert_almost_equal(csf2_1, 0.4, decimal=2)

    # csf2_2 = mu1.companion_star_fraction(70.0)
    # np.testing.assert_almost_equal(csf2_2, 2.0, decimal=2)

    # csf2_3 = mu1.companion_star_fraction(0.1)
    # np.testing.assert_almost_equal(csf2_3, 0.159, decimal=2)

    # Test brown dwarf csf
    csf_bd1 = mu1.companion_star_fraction(0.07)
    csf_bd2 = mu1.companion_star_fraction(0.04)
    csf_bd3 = mu1.companion_star_fraction(0.01)
    assert np.isclose(csf_bd1, 0.16, atol=0.01)
    assert np.isclose(csf_bd2, 0.08, atol=0.01)
    assert np.isclose(csf_bd3, 0.0, atol=1e-6)


def test_resolvedmult():
    """
    Test creating a MultiplicityResolvedDK object 
    and that the parameters it's populated with are correct.
    Updated to test for specific brown dwarf characteristics.
    """
    from spisea import synthetic, evolution, atmospheres, reddening, ifmr
    from spisea.imf import imf, multiplicity
    
    # Fetch isochrone
    logAge = 6.70 # Age in log(years)
    AKs = 1.0 # Ks filter extinction in mags
    dist = 4000 # distance in parsecs
    metallicity = 0 # metallicity in [M/H]
    atm_func = atmospheres.get_merged_atmosphere
    evo_merged = evolution.MergedPhillipsBaraffePisaEkstromParsec()
    redlaw = reddening.RedLawCardelli(3.1) # Rv = 3.1
    filt_list = ['nirc2,J', 'nirc2,Kp']
    
    startTime = time.time()

    iso_merged = synthetic.IsochronePhot(logAge, AKs, dist, metallicity=metallicity,
                                 evo_model=evo_merged, atm_func=atm_func,
                                 filters=filt_list, red_law=redlaw,
                                 mass_sampling=3)
    print('Constructed isochrone: %d seconds' % (time.time() - startTime))
    
    # Now we can make the cluster. 
    clust_mtot = 10**4.
    clust_multiplicity = multiplicity.MultiplicityResolvedDK()

    # Multiplicity is defined in the IMF object
    clust_imf_Mult = imf.Salpeter_Kirkpatrick_2024(multiplicity=clust_multiplicity)
    
    # Make clusters
    clust_Mult = synthetic.ResolvedCluster(iso_merged, clust_imf_Mult, clust_mtot)

    clust_Mult_ss = clust_Mult.star_systems
    
    print('Constructed cluster: %d seconds' % (time.time() - startTime))
    
    #check if columns were created
    assert 'log_a' in clust_Mult.companions.colnames
    assert 'e' in clust_Mult.companions.colnames
    assert 'i' in clust_Mult.companions.colnames
    assert 'Omega' in clust_Mult.companions.colnames
    assert 'omega' in clust_Mult.companions.colnames
    
    #check values are in correct range
    assert all(10**i<= 2000 and 10**i>= 0 for i in clust_Mult.companions['log_a']) #max separation is 2000 AU
    assert all(i<= 1 and i>= 0 for i in clust_Mult.companions['e'])
    assert all(i<= 180 and i>= 0 for i in clust_Mult.companions['i'])
    assert all(i<= 360 and i>= 0 for i in clust_Mult.companions['omega'])
    assert all(i<= 360 and i>= 0 for i in clust_Mult.companions['Omega'])
    
    #checks sign for inclination is being randomly genarated
    assert any(i > 90 for i in clust_Mult.companions['i']) and any(i < 90 for i in clust_Mult.companions['i'])
    
    #checks eccentricity follows f(e) = 2e pdf
    n, bins = np.histogram(clust_Mult.companions['e'], density = True)
    bin_centers = 0.5*(bins[1:] + bins[:-1])
    assert all(np.abs(i) < 0.3 for i in 2*bin_centers - n)
    
    #checks shape of inclination histogram is sin(i)
    n, bins = np.histogram(clust_Mult.companions['i'])
    bin_centers = 0.5*(bins[1:] + bins[:-1])
    assert all(np.abs(i) < 0.15 for i in n/max(n) - np.sin(np.pi*bin_centers/180))

    #checks for brown dwarf specific features
    bd_idx = np.where(clust_Mult.star_systems['mass'] < 0.08)[0]

    #check there is only one possible companion per BD
    assert all(clust_Mult.star_systems['N_companions'][bd_idx] <= 1), \
    "Brown dwarf primaries have >1 companion."
    
    comp_rows = []
    start = 0
    for ii, N in enumerate(clust_Mult.star_systems['N_companions']):
        if ii in bd_idx and N > 0:
            comp_rows.extend(range(start, start+N))
        start += N
    
    bd_companions = clust_Mult.companions[comp_rows]
    
    if len(bd_companions) > 30:  # only test if enough BD binaries
        mean_log_a = np.mean(bd_companions['log_a'])
        std_log_a = np.std(bd_companions['log_a'])
    
        #expect lognormal centered near log10(2.9 AU), width ~0.21
        assert abs(mean_log_a - np.log10(2.9)) < 0.25, \
            f"BD mean log(a) off: {mean_log_a:.2f}"
        assert abs(std_log_a - 0.21) < 0.15, \
            f"BD sigma log(a) off: {std_log_a:.2f}"
    
    return
