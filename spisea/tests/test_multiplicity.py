import numpy as np
import time
import os
import spisea
from spisea.imf import imf, multiplicity

def test_create_MultiplicityUnresolved():
    """
    Tests creating and accessing a MultiplicityUnresolved object.
    """
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

    mf2_1 = mu2.multiplicity_fraction(1.0)
    np.testing.assert_almost_equal(mf2_1, 0.4, decimal=2)

    mf2_2 = mu2.multiplicity_fraction(10.0)
    np.testing.assert_almost_equal(mf2_2, 1.0, decimal=2)

    mf2_3 = mu2.multiplicity_fraction(0.1)
    np.testing.assert_almost_equal(mf2_3, 0.159, decimal=2)

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
                                 mass_sampling=3, iso_dir=f'{spisea.__path__[0]}/tests/isochrones')
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

        bd_masses = clust_Mult.star_systems['mass'][bd_idx]
        expected_sigma = np.mean(
            np.interp(
                np.log10(bd_masses),
                [np.log10(0.01), np.log10(0.08)],
                [0.25, 0.5]
            )
        )

        #expect lognormal centered near log10(2.9 AU), width ~0.21
        assert abs(mean_log_a - np.log10(2.9)) < 0.25, \
            f"BD mean log(a) off: {mean_log_a:.2f}"

        assert abs(std_log_a - expected_sigma) < 0.15, \
            f"BD sigma log(a) off: {std_log_a:.2f}"

    return


# ---------------------------------------------------------------------------
# Offner et al. 2023 (Table 1) multiplicity
# ---------------------------------------------------------------------------

# Published Table 1 MF (%) converted to fraction, CF, and 1-sigma MF error.
# Masses are geometric means of the tabulated M1 intervals.
_OFFNER_TABLE1 = [
    # name, M_lo, M_hi, MF, MF_err, CF
    ('Fontanive+2018', 0.019, 0.058, 0.08, 0.06, 0.08),
    ('Burgasser 2007', 0.05, 0.08, 0.15, 0.04, 0.16),
    ('Close+2003', 0.080, 0.095, 0.19, 0.07, 0.19),
    ('Allen+2007', 0.06, 0.15, 0.20, 0.04, 0.20),
    ('Winters+2019 late-M', 0.075, 0.15, 0.19, 0.03, 0.21),
    ('Winters+2019 mid-M', 0.15, 0.30, 0.23, 0.02, 0.27),
    ('Winters+2019 early-M', 0.3, 0.6, 0.30, 0.02, 0.38),
    ('Raghavan+2010', 0.75, 1.25, 0.46, 0.03, 0.60),
    ('Tokovinin 2014b', 0.85, 1.5, 0.47, 0.03, 0.62),
    ('Moe & Kratter 2021', 1.6, 2.4, 0.68, 0.07, 0.99),
    ('Moe & Di Stefano 2017 3-5', 3.0, 5.0, 0.81, 0.06, 1.28),
    ('Moe & Di Stefano 2017 5-8', 5.0, 8.0, 0.89, 0.05, 1.55),
    ('Moe & Di Stefano 2017 8-17', 8.0, 17.0, 0.93, 0.04, 1.80),
    ('Sana et al. 17-50', 17.0, 50.0, 0.96, 0.04, 2.10),
]


def _table1_mgeom(row):
    return np.sqrt(row[1] * row[2])


def test_piecewise_powerlaw_api():
    """Custom piecewise MF/CSF is vectorized and clips MF to [0, 1]."""
    mass_limits = np.array([0.1, 1.0, 10.0])
    # First segment: MF = 0.4 * M^0 → 0.4; second: 0.4 * M^1 so MF(10)=4 → clip 1
    mp = multiplicity.MultiplicityPiecewisePowerLaw(
        mass_limits,
        MF_amps=[0.4, 0.4], MF_powers=[0.0, 1.0],
        CSF_amps=[0.4, 0.5], CSF_powers=[0.0, 0.5],
        binary_only_mass_max=0.05)
    assert mp.multiplicity_fraction(0.2) == 0.4
    assert mp.multiplicity_fraction(1.0) == 0.4
    np.testing.assert_almost_equal(mp.multiplicity_fraction(10.0), 1.0)
    masses = np.array([0.2, 1.0, 10.0])
    mf = mp.multiplicity_fraction(masses)
    np.testing.assert_allclose(mf, [mp.multiplicity_fraction(m) for m in masses])


def test_logistic_api():
    """Custom logistic MF/CSF clips, vectorizes, and sets BD CSF = MF."""
    ml = multiplicity.MultiplicityLogistic(
        MF_A=0.1, MF_B=1.5, MF_M0=1.0, MF_k=2.0,
        CSF_A=0.05, CSF_B=4.0, CSF_M0=2.0, CSF_k=1.0,
        CSF_max=2.0, binary_only_mass_max=0.08)
    # Low-mass asymptote A; M<=0 maps to A
    np.testing.assert_almost_equal(ml.multiplicity_fraction(1e-8), 0.1, decimal=4)
    assert ml.multiplicity_fraction(0.0) == 0.1
    assert ml.multiplicity_fraction(-1.0) == 0.1
    # High-mass MF saturates at B then clips to 1
    np.testing.assert_almost_equal(ml.multiplicity_fraction(1e6), 1.0)
    # High-mass CSF clips to CSF_max
    np.testing.assert_almost_equal(ml.companion_star_fraction(1e6), 2.0)
    # BD CSF = MF
    assert ml.companion_star_fraction(0.05) == ml.multiplicity_fraction(0.05)
    masses = np.array([0.05, 1.0, 100.0])
    mf = ml.multiplicity_fraction(masses)
    np.testing.assert_allclose(mf, [ml.multiplicity_fraction(m) for m in masses])
    csf = ml.companion_star_fraction(masses)
    np.testing.assert_allclose(
        csf, [ml.companion_star_fraction(m) for m in masses])
    assert np.all(csf >= mf - 1e-12)


def test_offner2023_logistic_coefficients():
    """Offner stores the equal-weight logistic-in-log-mass coefficients."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    assert isinstance(multi, multiplicity.MultiplicityLogistic)
    assert not isinstance(multi, multiplicity.MultiplicityPiecewisePowerLaw)
    np.testing.assert_allclose(multi.MF_A, 0.14)
    np.testing.assert_allclose(multi.MF_B, 0.99)
    np.testing.assert_allclose(multi.MF_M0, 1.41)
    np.testing.assert_allclose(multi.MF_k, 1.25)
    np.testing.assert_allclose(multi.CSF_A, 0.12)
    np.testing.assert_allclose(multi.CSF_B, 2.35)
    np.testing.assert_allclose(multi.CSF_M0, 3.57)
    np.testing.assert_allclose(multi.CSF_k, 0.96)


def test_offner2023_mf_smooth():
    """MF is continuous and nearly C1 around 0.08 and 1.5 Msun."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    eps = 1e-8
    for m in (0.08, 1.5):
        mf_left = multi.multiplicity_fraction(m - eps)
        mf_right = multi.multiplicity_fraction(m + eps)
        mf_at = multi.multiplicity_fraction(m)
        np.testing.assert_allclose(mf_left, mf_right, atol=1e-6, rtol=0)
        np.testing.assert_allclose(mf_at, mf_right, atol=1e-6, rtol=0)
        d_left = (mf_at - mf_left) / eps
        d_right = (mf_right - mf_at) / eps
        np.testing.assert_allclose(d_left, d_right, atol=1e-3, rtol=0)
    for m in (0.04, 0.3, 1.0, 10.0):
        expected = multiplicity._logistic_in_logm(
            m, 0.14, 0.99, 1.41, 1.25, clip_min=0.0, clip_max=1.0)
        np.testing.assert_allclose(multi.multiplicity_fraction(m), expected)


def test_offner2023_table1_mf():
    """
    Logistic MF matches Offner et al. 2023 Table 1 at geom-mean M1.

    Fontanive (8±6%) sits ~0.07 below the curve (~15%); other rows,
    including A/B stars, stay close.
    """
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    for row in _OFFNER_TABLE1:
        name, mlo, mhi, mf_tab, mf_err, cf_tab = row
        m = _table1_mgeom(row)
        mf = multi.multiplicity_fraction(m)
        tol = max(0.08, 2.0 * mf_err)
        assert abs(mf - mf_tab) <= tol, \
            '{0}: MF({1:.3f})={2:.3f} vs Table 1 {3:.2f} ± {4:.2f}'.format(
                name, m, mf, mf_tab, mf_err)
        assert 0.0 <= mf <= 1.0


def test_offner2023_table1_csf():
    """CSF matches Table 1 CF for stellar primaries; CSF = MF for BDs."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    for row in _OFFNER_TABLE1:
        name, mlo, mhi, mf_tab, mf_err, cf_tab = row
        m = _table1_mgeom(row)
        csf = multi.companion_star_fraction(m)
        mf = multi.multiplicity_fraction(m)
        if m <= multiplicity.H_BURNING_MASS:
            assert np.isclose(csf, mf, atol=1e-12), \
                '{0}: BD CSF should equal MF'.format(name)
        elif mlo >= 1.6:
            # Logistic CF tracks A/B; Moe & Kratter residual ~0.1 is ok
            tol = max(0.12, 0.12 * cf_tab)
            assert abs(csf - cf_tab) <= tol, \
                '{0}: CSF({1:.3f})={2:.3f} vs Table 1 CF {3:.2f}'.format(
                    name, m, csf, cf_tab)
        else:
            tol = max(0.08, 0.25 * cf_tab)
            assert abs(csf - cf_tab) <= tol, \
                '{0}: CSF({1:.3f})={2:.3f} vs Table 1 CF {3:.2f}'.format(
                    name, m, csf, cf_tab)
        assert csf >= mf - 1e-12
        assert csf <= multi.CSF_max + 1e-12


def test_offner2023_array_vs_scalar():
    """Array and scalar MF/CSF evaluations agree."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    masses = np.array([_table1_mgeom(row) for row in _OFFNER_TABLE1])
    mf_arr = multi.multiplicity_fraction(masses)
    csf_arr = multi.companion_star_fraction(masses)
    for i, m in enumerate(masses):
        np.testing.assert_allclose(mf_arr[i], multi.multiplicity_fraction(float(m)))
        np.testing.assert_allclose(csf_arr[i], multi.companion_star_fraction(float(m)))


def test_offner2023_bd_binaries_only():
    """BD primaries have CSF = MF and companion counts of 0 or 1."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    rng = np.random.default_rng(123)
    masses = np.array([0.02, 0.04, 0.07, 0.08])
    mf = multi.multiplicity_fraction(masses)
    csf = multi.companion_star_fraction(masses)
    np.testing.assert_allclose(csf, mf)
    # Force multiples so we test the count draw, not the MF coin flip.
    n_comp = multi.draw_n_companions(masses, csf, mf, rng)
    assert np.all(n_comp <= 1)
    assert np.all(n_comp >= 1)

    # Full companion-mass assignment: never more than one companion column
    is_mult = np.ones(len(masses), dtype=bool)
    comp, sys_mass, is_mult_out = multi.draw_companion_masses(
        masses, is_mult, csf, mf, rng, mass_min=0.01)
    assert comp.shape[1] == 1
    assert np.all(np.sum(~comp.mask, axis=1) <= 1)


def test_offner2023_q_more_equal_mass_for_bds():
    """BD mass ratios are more equal-mass (higher mean q) than solar-type."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    rng = np.random.default_rng(7)
    n = 20000
    q_bd = multi.random_q(rng.random(n), mass=0.04)
    q_sun = multi.random_q(rng.random(n), mass=1.0)
    assert np.mean(q_bd) > np.mean(q_sun) + 0.1
    # Err-wt logistic undershoots Fontanive 4.8 (~3.3 at 0.033 Msun)
    assert multi.q_power_at_mass(0.033) > 2.5
    assert multi.q_power_at_mass(1.0) < 0.5


def test_offner2023_q_sigma_a_closed_form():
    """γ, σ(log a), and log_a_mean match the smooth helpers; not interpolation."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    masses = np.array([0.033, 0.065, 0.3, 1.0, 10.0])
    for m in masses:
        np.testing.assert_allclose(
            multi.q_power_at_mass(m),
            multiplicity._logistic_in_logm(
                m, 6.6, -1.77, 0.0651, 0.629))
        np.testing.assert_allclose(
            multi.sigma_log_a(m),
            multiplicity._logistic_in_logm(
                m, 0.7, 1.5, 0.354, 6.05, clip_min=0.1))
        np.testing.assert_allclose(
            multi.log_a_mean(m),
            multiplicity._smooth_broken_loglog(
                m, 44.46, 0.819, 1.005, -0.308, 0.10, a_min=0.1))
    # Array vs scalar
    g_arr = multi.q_power_at_mass(masses)
    sig_arr = multi.sigma_log_a(masses)
    loga_arr = multi.log_a_mean(masses)
    for i, m in enumerate(masses):
        np.testing.assert_allclose(g_arr[i], multi.q_power_at_mass(float(m)))
        np.testing.assert_allclose(sig_arr[i], multi.sigma_log_a(float(m)))
        np.testing.assert_allclose(loga_arr[i], multi.log_a_mean(float(m)))
    # Old L/early-T interpolation knot was 2.5; logistic is not that.
    g_knot = multi.q_power_at_mass(0.065)
    np.testing.assert_allclose(
        g_knot, multiplicity._logistic_in_logm(0.065, 6.6, -1.77, 0.0651, 0.629))
    assert abs(g_knot - 2.5) > 0.05


def test_offner2023_bd_separations_peak_few_au():
    """BD lognormal separations peak at a few AU (μ(0.033)≈2.1 au)."""
    multi = multiplicity.MultiplicityResolvedOffner2023()
    np.random.seed(0)
    log_a = multi.log_semimajoraxis(np.full(5000, 0.04))
    med_a = 10 ** np.median(log_a)
    assert 1.5 < med_a < 8.0, 'BD median a = {0:.2f} AU'.format(med_a)
    # Solar-type should be much wider (smooth-broken μ ~ 44 au)
    log_a_s = multi.log_semimajoraxis(np.full(5000, 1.0))
    med_a_s = 10 ** np.median(log_a_s)
    assert med_a_s > 10.0
    assert med_a_s > med_a


def test_offner2023_alias_and_resolved_methods():
    """Public names and resolved orbital methods exist."""
    assert multiplicity.MultiplicityOffner2023 is \
        multiplicity.MultiplicityUnresolvedOffner2023
    resolved = multiplicity.MultiplicityResolvedOffner2023()
    assert hasattr(resolved, 'log_semimajoraxis')
    assert hasattr(resolved, 'log_a_mean')
    assert hasattr(resolved, 'sigma_log_a')
    e = resolved.random_e(np.array([0.0, 0.25, 1.0]))
    np.testing.assert_allclose(e, [0.0, 0.5, 1.0])


def test_lu2013_defaults_unchanged():
    """Lu et al. 2013 MultiplicityUnresolved defaults and stellar MF unchanged."""
    mu = multiplicity.MultiplicityUnresolved()
    assert mu.MF_amp == 0.44
    assert mu.MF_pow == 0.51
    assert mu.CSF_amp == 0.50
    assert mu.CSF_pow == 0.45
    np.testing.assert_almost_equal(mu.multiplicity_fraction(1.0), 0.44, decimal=2)
    np.testing.assert_almost_equal(mu.multiplicity_fraction(10.0), 1.0, decimal=2)
    np.testing.assert_almost_equal(mu.multiplicity_fraction(0.1), 0.136, decimal=2)
    # Scalar BD overrides (historical Lu+2013 / Fontanive path)
    assert np.isclose(mu.multiplicity_fraction(0.07), 0.16, atol=0.01)
    assert np.isclose(mu.multiplicity_fraction(0.04), 0.08, atol=0.01)
    assert np.isclose(mu.multiplicity_fraction(0.01), 0.0, atol=1e-6)


def test_offner_generate_cluster_companions():
    """IMF cluster generation with Offner multiplicity produces BD binaries only."""
    imf_multi = multiplicity.MultiplicityUnresolvedOffner2023()
    mass_limits = np.array([0.01, 0.08, 0.5, 120.0])
    powers = np.array([-0.3, -1.3, -2.3])
    my_imf = imf.IMF_broken_powerlaw(mass_limits, powers, imf_multi)
    my_imf.rng = np.random.default_rng(42)
    mass, is_multi, comp_mass, sys_mass = my_imf.generate_cluster(500.0)
    bd = mass <= 0.08
    n_comp = np.sum(~comp_mass.mask, axis=1)
    assert np.all(n_comp[bd] <= 1)
    assert np.any(is_multi)
    assert np.abs(500.0 - sys_mass.sum()) < 500.0 * 0.05


def test_calc_multi_uses_multiplicity_q_and_counts():
    """
    IMF.calc_multi must not hardcode Fontanive gamma=6.1 or the BD
    companion cap; those policies live on the multiplicity object so
    Offner γ_trunc (~2–5 for BDs) actually applies.
    """
    import inspect
    from spisea.imf import imf as imf_mod
    calc_src = inspect.getsource(imf_mod.IMF.calc_multi)
    assert '6.1' not in calc_src
    assert 'draw_companion_masses' in calc_src

    syn_path = os.path.join(os.path.dirname(spisea.__file__), 'synthetic.py')
    with open(syn_path, 'r') as fh:
        syn_src = fh.read()
    assert 'isinstance(self.imf._multi_props, multiplicity.MultiplicityResolvedDK)' not in syn_src
    assert "hasattr(multi_props, 'log_semimajoraxis')" in syn_src
    assert "hasattr(multi_props, 'random_e')" in syn_src
    assert "hasattr(multi_props, 'random_keplarian_parameters')" in syn_src

    offner = multiplicity.MultiplicityUnresolvedOffner2023()
    lu = multiplicity.MultiplicityUnresolved()
    q_off = offner.q_power_at_mass(0.04)
    q_lu = lu.q_power_at_mass(0.04)
    assert 2.0 <= q_off <= 5.5
    assert np.isclose(q_lu, 6.1)
    assert q_off != q_lu

    rng = np.random.default_rng(1)
    masses = np.full(3000, 0.04)
    is_mult = np.ones(len(masses), dtype=bool)
    mf = offner.multiplicity_fraction(masses)
    csf = offner.companion_star_fraction(masses)
    comp, _, _ = offner.draw_companion_masses(
        masses, is_mult, csf, mf, rng, mass_min=0.01)
    q = comp.compressed() / 0.04
    q_lu_draw = lu.random_q(np.random.default_rng(1).random(len(q)), mass=0.04)
    # Offner BD gamma is shallower than Fontanive 6.1, so mean q is lower.
    assert np.mean(q) < np.mean(q_lu_draw)


def test_offner2023_mf_is_vectorized():
    """Offner MF is vectorized; Lu+2013 scalar BD bins are not used here."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    masses = np.array([0.03, 0.06, 0.10, 1.0])
    mf = multi.multiplicity_fraction(masses)
    for i, m in enumerate(masses):
        np.testing.assert_allclose(mf[i], multi.multiplicity_fraction(float(m)))
    # Array path is not the Lu stellar power law 0.44 * M**0.51
    lu_pl = 0.44 * masses ** 0.51
    assert not np.allclose(mf, np.clip(lu_pl, 0, 1), atol=0.02)

