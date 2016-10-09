import numpy as np
import nose.tools

def test_prim_power():
    from .. import imf

    mass_limits = np.array([0.1, 1.0, 100.0])
    powers = np.array([-1.8, -2.0])

    imf_tmp = imf.IMF_broken_powerlaw(mass_limits, powers)

    (mass, is_multi, c_mass, s_mass) = imf_tmp.generate_cluster(1e3)

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
    print('Mass Sum: ', mass.sum(), '  (should be 1e3)')

    return
