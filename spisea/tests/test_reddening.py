import numpy as np
from spisea import reddening, synthetic
import pdb


def test_RedLawBrokenPowerLaw(plots=False):
    # Start with 2-segment law
    lambda_limits = [2.14, 1.63, 1.27]
    alpha_vals = [2.23, 2.44]
    K_wave = 2.14

    red_law = reddening.RedLawBrokenPowerLaw(lambda_limits, alpha_vals, K_wave)

    # Test redlaw at several wavelengths to amek sure it is what we expect
    
    # If desired (debug only), make plot to see what law looks like
    if plots:

        pdb.set_trace()


    

    pdb.set_trace()
    return

