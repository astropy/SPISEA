import numpy as np
from spisea import reddening, synthetic
import pdb


def test_RedLawBrokenPowerLaw(plots=False):
    # Start with 2-segment law
    lambda_limits = [2.16, 1.63, 1.27]
    alpha_vals = [2.23, 2.44]
    K_wave = 2.14

    red_law = reddening.RedLawBrokenPowerLaw(lambda_limits, alpha_vals, K_wave)

    # Calculate what the redlaw *should* be across wavelength range,
    # by-hand
    wave_test = np.arange(1.27, 2.15+0.001, 0.01)
    law_test = np.ones(len(wave_test)) * np.nan

    # 2.14 - 1.63
    idx = np.where( (wave_test > 1.63) & (wave_test <= 2.17))
    coeff = 1
    law_test[idx] = coeff * wave_test[idx] ** (-2.23)

    # 1.63 - 1.27
    idx = np.where( (wave_test >= 1.27) & (wave_test <= 1.63))
    coeff  = (1.63 ** -2.23) / (1.63 ** 2.44)
    law_test[idx] = coeff * wave_test[idx] ** (-2.23)
    assert np.sum(np.isnan(law_test)) == 0

    # Put in terms of A_lambda / A_Ks, like the reddening object
    idx_K = np.where(abs(wave_test - K_wave) == min(abs(wave_test - K_wave)))
    law_test /= law_test[idx_K]

    # Compare law_test and the output from the redlaw object
    law_output = red_law.broken_powerlaw(wave_test, 1)
    
    assert len(law_test) == len(law_output)
    assert np.sum(np.isnan(law_output)) == 0

    # Calculate the difference between test calc and code output
    diff = law_output - law_test
    
    # If desired (debug only), make plot to see what law looks like
    if plots:

        pdb.set_trace()


    

    pdb.set_trace()
    return

