import numpy as np
import nose.tools
    
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


    
