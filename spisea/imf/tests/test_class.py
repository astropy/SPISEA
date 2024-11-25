import numpy as np
import time
import pdb

def test_Muzic():
    from .. import imf
    from .. import multiplicity
    
    # Make multiplicity object
    imf_multi = multiplicity.MultiplicityUnresolved()

    # Make IMF object; we'll use a broken power law with the parameters from Muzic+17
    my_imf = imf.Muzic_2017(multiplicity=imf_multi)

    # Define total cluster mass
    M_cl = 10**5.

    mass, isMulti, compMass, sysMass = my_imf.generate_cluster(M_cl)

    # Make sure that the total mass is always within the expected
    # range of the requested mass.
    assert np.abs(M_cl - sysMass.sum()) < 120.0

    # Check that enough companions were generated.
    # Should be greater than 25% of the stars with companions.
    mdx = np.where(isMulti)[0]
    for mm in mdx:
        assert len(compMass[mm]) > 0

    return

def test_CombinedWeidnerKroupaKirkpatrick_2024():
    from .. import imf
    from .. import multiplicity
    
    # Make multiplicity object
    imf_multi = multiplicity.MultiplicityUnresolved()

    # Make IMF object; we'll use a broken power law with the parameters from Muzic+17
    my_imf = imf.CombinedWeidnerKroupaKirkpatrick_2024(multiplicity=imf_multi)

    # Define total cluster mass
    M_cl = 10**5.

    mass, isMulti, compMass, sysMass = my_imf.generate_cluster(M_cl)

    # Make sure that the total mass is always within the expected
    # range of the requested mass.
    assert np.abs(M_cl - sysMass.sum()) < 120.0

    # Check that enough companions were generated.
    # Should be greater than 25% of the stars with companions.
    mdx = np.where(isMulti)[0]
    for mm in mdx:
        assert len(compMass[mm]) > 0

    return