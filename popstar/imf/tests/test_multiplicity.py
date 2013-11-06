def test_create_MultiplicityUnresolved():
    """
    Tests creating and accessing a MultiplicityUnresolved object.
    """
    from .. import multiplicity

    # All default parameters -- check their values
    mu1 = multiplicity.MultiplicityUnresolved()
    assert mu1.MF_amp == 0.44
    assert mu1.MF_power == 0.51
    assert mu1.CSF_amp == 0.50
    assert mu1.CSF_power == 0.45
    assert mu1.CSF_max == 3
    assert mu1.q_power == -0.4
    assert mu1.q_min == 0.01

    # Test setting different parameters
    mu2 = multplicity.MultiplicityUnresolved(MF_amp=0.4, 
                                             MF_power=0.4,
                                             CSF_amp=0.4, 
                                             CSF_power=0.4, 
                                             CSF_max=4,
                                             q_power=0.4, 
                                             q_min=0.04)
    assert mu2.MF_amp == 0.4 
    assert mu2.MF_power == 0.4
    assert mu2.CSF_amp == 0.4 
    assert mu2.CSF_power == 0.4 
    assert mu2.CSF_max == 4
    assert mu2.q_power == 0.4 
    assert mu2.q_min == 0.04

    

    
