import numpy as np
import nose.tools
import time
import pdb

def test_generate_cluster():
    from .. import imf
    from .. import multiplicity
    
    # Make multiplicity object
    imf_multi = multiplicity.MultiplicityUnresolved()

    # Make IMF object; we'll use a broken power law with the parameters from Kroupa+01
    massLimits = np.array([0.08, 0.5, 1, 120]) # Define boundaries of each mass segement
    powers = np.array([-1.3, -2.3, -2.3]) # Power law slope associated with each mass segment
    my_imf = imf.IMF_broken_powerlaw(massLimits, powers, imf_multi)

    # Define total cluster mass
    M_cl = 10**5.

    mass, isMulti, compMass, sysMass = my_imf.generate_cluster(M_cl)

    # Make sure that the total mass is always within the expected
    # range of the requested mass.
    assert np.abs(M_cl - sysMass.sum()) < 120.0

    return


def test_prim_power():
    from .. import imf

    #mass_limits = np.array([0.1, 1.0, 100.0])
    #powers = np.array([-2.0, -1.8])
    mass_limits = np.array([1.0, 100.0])
    powers = np.array([-2.0])

    imf_tmp = imf.IMF_broken_powerlaw(mass_limits, powers)

    Mcl = 1e5
    (mass, is_multi, c_mass, s_mass) = imf_tmp.generate_cluster(Mcl)

    print('mass shape = ', mass.shape)
    print('mass array:')
    print(mass)
    print('is_multi array:')
    print(is_multi)
    print('c_mass array:')
    print(c_mass)
    print('s_mass array:')
    print(s_mass)
    print('np.isfinite(mass):')
    print(np.isfinite(mass))
    print('Mass Sum: ', mass.sum(), '  (should be {0})'.format(Mcl))

    return

def test_xi():
    from .. import imf

    import cProfile, pstats, io
    #from pstats import SortKey
    pr = cProfile.Profile()

    mass_limits = np.array([0.1, 1.0, 10.0, 100.0])
    powers = np.array([-0.3, -1.5, -2.3])
    imf_tmp = imf.IMF_broken_powerlaw(mass_limits, powers)

    ##########
    #
    # Test validity of returned values. 
    # 
    ##########
    N_size = 10
    m = np.linspace(0.2, 20, N_size)
    val_good = np.array([1.6206566 , 0.26895718, 0.10135922, 0.05639448, 0.03703704,
                         0.0243668 , 0.01613091, 0.01137139, 0.00839526, 0.00642142])
    
    val_test = np.zeros(len(m), dtype=float)
    for ii in range(len(m)):
        val_test[ii] = imf_tmp.xi(m[ii])

    # print(repr(val_test))
    np.testing.assert_almost_equal(val_test, val_good)

    ##########
    #
    # Performance testing
    # 
    ##########
    t1 = time.time()
    # pr.enable()
    
    # Run a time test
    N_size = int(1e4)
    m = np.random.uniform(1.1, 99.0, size=N_size)
    foo1 = imf_tmp.xi(m)
    
    # pr.disable()
    t2 = time.time()

    print('test_xi() runtime = {0:.3f} s for {1:d} masses'.format(t2 - t1, N_size))
    
    # s = io.StringIO()
    # sortby = SortKey.CUMULATIVE
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print(s.getvalue())    

    return

def test_mxi():
    from .. import imf

    import cProfile, pstats, io
    #from pstats import SortKey
    pr = cProfile.Profile()

    mass_limits = np.array([0.1, 1.0, 10.0, 100.0])
    powers = np.array([-0.3, -1.5, -2.3])
    imf_tmp = imf.IMF_broken_powerlaw(mass_limits, powers)

    ##########
    #
    # Test validity of returned values. 
    # 
    ##########
    N_size = 10
    m = np.linspace(0.2, 20, N_size)
    val_good = np.array([0.32413132, 0.64549722, 0.4662524 , 0.38348249, 0.33333333,
                         0.27290819, 0.21615424, 0.17739363, 0.14943557, 0.12842838])
    val_test = np.zeros(len(m), dtype=float)
    for ii in range(len(m)):
        val_test[ii] = imf_tmp.m_xi(m[ii])

    # print(repr(val_test))
    np.testing.assert_almost_equal(val_test, val_good)
    assert val_test.shape == m.shape

    ##########
    #
    # Performance testing
    # 
    ##########
    t1 = time.time()
    # pr.enable()
    
    # Run a time test
    N_size = int(1e4)
    m = np.random.uniform(1.1, 99.0, size=N_size)
    foo1 = imf_tmp.m_xi(m)
    assert foo1.shape == m.shape

    # pr.disable()
    t2 = time.time()

    print('test_mxi() runtime = {0:.3f} s for {1:d} masses'.format(t2 - t1, N_size))
    
    # s = io.StringIO()
    # sortby = SortKey.CUMULATIVE
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print(s.getvalue())    

    return

def test_theta_closed():
    from .. import imf

    import cProfile, pstats, io
    #from pstats import SortKey
    
    mass_limits = np.array([0.1, 1.0, 10.0, 100.0])
    powers = np.array([-0.3, -1.5, -2.3])
    imf_tmp = imf.IMF_broken_powerlaw(mass_limits, powers)

    ##########
    #
    # Test validity of returned values. 
    # 
    ##########
    N_size = 10
    m = np.linspace(0.2, 20, N_size)
    val_good = np.array([[1., 0., 0.],
       [1., 1., 0.],
       [1., 1., 0.],
       [1., 1., 0.],
       [1., 1., 0.],
       [1., 1., 1.],
       [1., 1., 1.],
       [1., 1., 1.],
       [1., 1., 1.],
       [1., 1., 1.]])

    val_test = np.zeros((len(m), len(powers)), dtype=float)
    for ii in range(len(m)):
        val_test[ii] = imf.theta_closed(m[ii] - imf_tmp._m_limits_low)

    np.testing.assert_equal(val_test, val_good)

    

    ##########
    #
    # Speed tests and performance profiling.
    # 
    ##########
    N_size = 10000
    m = np.linspace(1.1, 99, N_size)
    
    tmp = np.zeros((len(m), len(powers)), dtype=float)

    t1 = time.time()
    # pr = cProfile.Profile()
    # pr.enable()
    
    for ii in range(len(m)):
        tmp[ii] = imf.theta_closed(m[ii] - imf_tmp._m_limits_low)

    # pr.disable()
    t2 = time.time()

    print('Runtime = {0:.3f} s for {1:d} masses'.format(t2 - t1, N_size))
    
    # s = io.StringIO()
    # sortby = SortKey.CUMULATIVE
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print(s.getvalue())


    return
    
