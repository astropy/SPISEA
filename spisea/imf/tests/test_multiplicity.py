import numpy as np
import nose.tools
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
    nose.tools.assert_almost_equal(mf1_1, 0.44, places=2)

    mf1_2 = mu1.multiplicity_fraction(10.0)
    nose.tools.assert_almost_equal(mf1_2, 1.0, places=2)

    mf1_3 = mu1.multiplicity_fraction(0.1)
    nose.tools.assert_almost_equal(mf1_3, 0.136, places=2)

    # Second set of multiplicity parameters
    mu2 = multiplicity.MultiplicityUnresolved(MF_amp=0.4, MF_power=0.4,
                                             CSF_amp=0.4, CSF_power=0.4, CSF_max=4,
                                             q_power=0.4, q_min=0.04)

    mf2_1 = mu1.multiplicity_fraction(1.0)
    nose.tools.assert_almost_equal(mf2_1, 0.44, places=2)

    mf2_2 = mu1.multiplicity_fraction(10.0)
    nose.tools.assert_almost_equal(mf2_2, 1.0, places=2)

    mf2_3 = mu1.multiplicity_fraction(0.1)
    nose.tools.assert_almost_equal(mf2_3, 0.136, places=2)


def test_multiplicity_fraction_array():
    """
    Test multiplicity_fraction() on the MultiplicityUnresolved object
    where the inputs and outputs are in array form.
    """
    from spisea.imf import multiplicity
    
    # First set of multiplicity parameters
    mu1 = multiplicity.MultiplicityUnresolved()

    mass_array = np.array([1.0, 10.0, 0.1])
    mf_array = mu1.multiplicity_fraction(mass_array)

    nose.tools.assert_almost_equal(mf_array[0], 0.44, places=2)
    nose.tools.assert_almost_equal(mf_array[1], 1.0, places=2)
    nose.tools.assert_almost_equal(mf_array[2], 0.136, places=2)
    
    
def test_companion_star_fraction():
    """
    Test the companion_star fraction on the MultiplicityUnresolved object.
    """
    from spisea.imf import multiplicity

    # First set of multiplicity parameters
    mu1 = multiplicity.MultiplicityUnresolved()

    csf1_1 = mu1.companion_star_fraction(1.0)
    nose.tools.assert_almost_equal(csf1_1, 0.5, places=2)

    csf1_2 = mu1.companion_star_fraction(70.0)
    nose.tools.assert_almost_equal(csf1_2, 3.0, places=2)

    csf1_3 = mu1.companion_star_fraction(0.1)
    nose.tools.assert_almost_equal(csf1_3, 0.177, places=2)

    # Second set of multiplicity parameters
    mu2 = multiplicity.MultiplicityUnresolved(MF_amp=0.4, MF_power=0.4,
                                              CSF_amp=0.4, CSF_power=0.4, CSF_max=2,
                                              q_power=0.4, q_min=0.04)

    # csf2_1 = mu1.companion_star_fraction(1.0)
    # nose.tools.assert_almost_equal(csf2_1, 0.4, places=2)

    # csf2_2 = mu1.companion_star_fraction(70.0)
    # nose.tools.assert_almost_equal(csf2_2, 2.0, places=2)

    # csf2_3 = mu1.companion_star_fraction(0.1)
    # nose.tools.assert_almost_equal(csf2_3, 0.159, places=2)


def test_resolvedmult():
    """
    Test creating a MultiplicityResolvedDK object 
    and that the parameters it's populated with are correct.
    """
    from spisea import synthetic, evolution, atmospheres, reddening, ifmr
    from spisea.imf import imf, multiplicity
    
    # Fetch isochrone
    logAge = 6.70 # Age in log(years)
    AKs = 1.0 # Ks filter extinction in mags
    dist = 4000 # distance in parsecs
    metallicity = 0 # metallicity in [M/H]
    atm_func = atmospheres.get_merged_atmosphere
    evo_merged = evolution.MISTv1()
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
    clust_imf_Mult = imf.Kroupa_2001(multiplicity=clust_multiplicity)
    
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
    
    return

    
