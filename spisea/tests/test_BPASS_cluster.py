import pytest
from spisea import synthetic
import loading_jlu_file
from loading_jlu_file import load_klf_by_radius
from spisea import imf
from spisea import ifmr
import numpy as np
from pytest import approx
import math
import matplotlib.pyplot as plt
py = plt

def test_third_figure_tests():
    """ In this test I will try to make sure that
    BPASS clusters' KLF's do not deviate too widely
    from that of the Observed KLF from the Young Galactic Center Stars
    from Lu et al. 2013. (around 6 million years of age, solar metallicity,
    8000 parsecs from the Earth, and 2.7 AKs extinction).
    For our purposes, we want to make sure that the Observed
    and BPASS counterpart's KLF's are no farther apart than 0.25.
    We may also inspect the KLF of the BPASS and MIST v.1. cluster
    and also the overall initial mass of both BPASS and MIST v1 counterparts
    to the MIST cluster """

    result = load_klf_by_radius(mask_for_log=True)
    magBin = result.Kp[1] - result.Kp[0]
    idx = np.where(result.Kp < 16)[0]
    klf_mag_bins = np.arange(9.0, 17, 1.0)
    binsKp = klf_mag_bins
    binEdges = binsKp[0:-1] + (binsKp[1:] - binsKp[0:-1]) / 2.0


    BPASS_iso = synthetic.Isochrone_Binary(6.78, 2.7,
                                           8000, 0.0,
                                           filters=['nirc2,Kp'])
   # If we want to use control group, let's use multiplicity = None.
    custom_IMF = imf.imf.IMF_broken_powerlaw(np.array([1, 150]),
                                             np.array([-1.7]),
                                             multiplicity=
                                             (imf.multiplicity.
                                              MultiplicityResolvedDK()))
    BPASS_Cluster = synthetic.Cluster_w_Binaries(BPASS_iso, custom_IMF,
                                                 1000000,
                                                 ifmr=ifmr.IFMR_Spera15())
    
    area = 116.098  # arcsec^2
    # Think of how many stars per area would there be in a cluster
    # similar to BPASS_Cluster
    # but with lower mass. We use scaling.
    scaleFactorBPS = ((17000 /
                       BPASS_Cluster.star_systems['systemMass'].sum()) /
                      area)
    MIST_iso = synthetic.IsochronePhot(6.78, 2.7, 8000, 0.0,
                                       filters=['nirc2,Kp'])
    MIST_Cluster = synthetic.ResolvedCluster(MIST_iso, custom_IMF,
                                             1000000,
                                         ifmr=ifmr.IFMR_Spera15())
    totl_mist =(MIST_Cluster.star_systems['m_nirc2_Kp']
            [np.where((MIST_Cluster.star_systems['isWR'] == 0))[0]])
    totl_BPASS = (BPASS_Cluster.star_systems['m_nirc2_Kp']
               [np.where((~BPASS_Cluster.star_systems['isWR']))[0]])
    mist_scale = ((17000 /
                   MIST_Cluster.star_systems['systemMass'].sum()) /
                  area)
    weightsMST = np.array([1.0 for x in totl_mist])
    weightsMST *= mist_scale
    weightsBPS = np.array([1.0 for x in totl_BPASS])
    weightsBPS *= scaleFactorBPS
    # Binning the K' magnitudes into their respective bins
    # And plotting KLF's of BPASS_Cluster, MIST clusters, and
    # Observed KLF
    # n, bins, patches are exactly the same objects that are
    # described as outputs of matplotlib.hist
    # Right below: BPASS cluster KLF histogram's bin values, bin-edges, and patches
    n, bins, patches = py.hist(totl_mist, bins=binEdges, histtype='step',
                               weights=weightsMST, color='green', label='MISTv.1 Model',
                               align='mid', linewidth=1.5)
    # Right below: BPASS cluster KLF histogram's bin values, bin-edges, and patches
    n2, bins2, patches2 = py.hist(totl_BPASS, bins=binEdges, histtype='step',
                                  weights=weightsBPS, color='blue',
                                  label='(BPASS Model KLF)',
                                  align='mid', linewidth=1.5)
    py.errorbar(result.Kp[idx], result.KLF_ext_cmp_sp_im_noWR[idx],
                fmt='ro-', xerr=magBin/2.0, capsize=0, linewidth=2)
    py.errorbar(result.Kp[idx] ,result.KLF_ext_cmp_sp_im_noWR[idx],
                fmt='ro-', yerr=result.eKLF_ext_cmp_sp_im_noWR[idx],
                linewidth=2,
                label='Observed')
    py.xlim(8.5, 15.5)
    py.xlabel('Kp magnitude')
    py.ylabel("stars / (arcsecond^2 mag)")
    py.title("KLF's at Age = %d Myr" % (10**(6.78 - 6)), fontsize=14)
    py.legend(loc='upper left', numpoints=1)
    py.savefig('Comparisons_w_Real_Data.png')
    # Now compare differences in the heights of BPASS cluster KLF and given KLF
    for x in range(len(binEdges)):
        assert np.abs(result.KLF_ext_cmp_sp_im_noWR[idx][x] - n[x]) == approx(0, abs=0.25)
    print("BPASS IMF does not seem too far off when compared to the BPASS IMF")
    # Now compare differences in the BPASS and MIST KLF's (or values
    # of bins in histograms for both KLF's)
    for x in range(len(binEdges)):
        assert np.abs(n2[x] - n[x]) == approx(0, abs=0.25)
    # Now compare the initial cluster mass of both the BPASS
    # and MIST cluster
    assert (np.abs(np.log10(BPASS_Cluster()) -
                   np.log10(MIST_Cluster())) == approx(0, rel=0.2))
    

def test_second_figure_BPASS():
    """ In this test, I will generate figures comparing the
    current mass distributions of noncompact remnant stars
    (non-phase = -99 stars) in both BPASS and MIST clusters
    of 10^9.2 years of age and 1/10th of solar metallicity.
    In the tests I will create and save a figure juxtaposing
    the current mass distributions.
    
    Things I will be testing:
    - Will there be no main sequence stars of 
    > ~ 5 solar masses (thanks Dr. McKee)?
    - Will the highest mass stars have 
    
    
    """
    # Check if the evolution class works fine
    # BPASS Isochrone used to create cluster
    iso1 = synthetic.Isochrone_Binary(9.2, 0.5, 2000,
                                  math.log10(0.1),
                                  mass_sampling=1)
    # MIST v1 Isochrone used to create cluster
    iso2 = synthetic.IsochronePhot(9.2, 0.5, 2000,
                                   math.log10(0.1),
                                   recomp=True)
    clus_1 = synthetic.Cluster_w_Binaries(iso1,
                                      imf.imf.IMFSalpeter1955(multiplicity=
                                                          imf.multiplicity.
                                                          MultiplicityResolvedDK()),
                                      200000, ifmr=ifmr.IFMR_Spera15())
    clus_2 = synthetic.ResolvedCluster(iso2,
                                       imf.imf.IMFSalpeter1955(multiplicity=
                                                           imf.multiplicity.
                                                           MultiplicityResolvedDK()),
                                       200000, ifmr=ifmr.IFMR_Spera15())
    # extract the BPASS isochrone's Primary stars' and companion stars' current mass
    star_systems = clus_1.star_systems
    companions = clus_1.companions
    prims_cm1 = (star_systems['mass_current']
                 [np.where(star_systems['phase'] == 5)[0]])
    companions_cm1 = (companions['mass_current']
                      [np.where(companions['phase'] == 5)[0]])
    plt.hist(np.hstack((prims_cm1, companions_cm1)), 40, (0, 10),
             density=True, histtype='step', label="BPASS")
    star_systems2 = clus_2.star_systems
    companions2 = clus_2.companions
    prims_cm2 = star_systems2['mass_current'][np.where(star_systems2['phase'] <
                                                   101)[0]]
    companions_cm2 = companions2['mass_current'][np.where(companions2['phase'] <
                                                      101)[0]]
    plt.hist(np.hstack((prims_cm2 , companions_cm2)), 40, (0, 10),
         density=True, histtype='step', label="MISTv1")
    plt.xlabel("Current mass of the star in solar masses +"
               " MISTv1 clusters (PDF) (binsize = 0.25 solar mass)")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    py.savefig('Comparing_BPASS_Current_Mass_Dist_w_MIST.png')
    
def test_figure_one_BPASS():
    """ In this test, I will be running the plotting code comparing BPASS
    and SPISEA codes. The BPASS and SPISEA isochrones are 10^8.2 years old and
    have solar metallicity. I will also check if the supposed helium stars seem
    to fit the profile of what is associated with that stellar type. In other words,
    I verify whether the stars that are in what I think is the helium star area
    (4.25<log Teff<4.75 and 27.5 < log L < 29) have initial mass of greater than
    the solar mass.
    """
    iso_bpass = synthetic.Isochrone_Binary(8.2, 0.0, 2000, 0.0)
    # New MIST v.1 isochrone for same metallicity
    iso2 = synthetic.IsochronePhot(8.2, 0.0, 2000, math.log(1), recomp=True)
    # Now I will generate the plot 
    plt.plot(np.log10(iso_bpass.primaries['Teff']),
             np.log10(iso_bpass.primaries["L"]), "r.")
    plt.plot(np.log10(iso_bpass.secondaries['Teff']),
             np.log10(iso_bpass.secondaries["L"]), "r.")
    plt.plot(np.log10(iso_bpass.singles['Teff']),
             np.log10(iso_bpass.singles["L"]), "r.", label="BPASS isochrone")
    plt.plot(np.log10(iso2.points['Teff']), np.log10(iso2.points["L"]), "b+", label="MISTv1", alpha=0.25)
    plt.xlabel("log10(T in kelvin)")
    plt.ylabel("log10(L in watts)")
    plt.title("HR Diagram of Isochrones at solar metallicity and 10**8.2 years age")
    plt.gca().invert_xaxis()
    plt.legend()
    plt.savefig("Comparing_BPASS_and_MIST_isochrones_8.2 log years.png")
    # Now as a sanity check, I will show that the possible helium stars in my
    # isochrone are more luminous, more hotter than the sun and are of higher
    # initial mass (hence are closer to the O-B type that we expect for these
    # types of stars).
    # Check such thing holds for primary stars
    He_star_region = iso_bpass.primaries[np.where((np.log10(iso_bpass.primaries['Teff']) <= 4.75) &
                                            (np.log10(iso_bpass.primaries['Teff']) >= 4.25) &
                                            (np.log10(iso_bpass.primaries['L']) >= 27.5) &
                                            (np.log10(iso_bpass.primaries['L']) <= 29))]
    # Do all of the primary stars in that region have initial mass
    # >= 2.0 solar mass 
    for_all = all(He_star_region['mass'] >= 2.0)
    assert for_all
    # see if any of the secondary stars have become helium stars
    He_star_second = iso_bpass.secondaries[np.where((np.log10(iso_bpass.primaries['Teff']) <= 4.75) &
                                                    (np.log10(iso_bpass.primaries['Teff']) >= 4.25) &
                                                    (np.log10(iso_bpass.primaries['L']) >= 27.5) &
                                                    (np.log10(iso_bpass.primaries['L']) <= 29))]
    for_all_2 = all(He_star_second['mass'] >= 2.0)
    assert for_all_2
    # To do: write subtest for checking the main sequence of the isochrone

def test_isochrones_dim_stars():
    """ In this test, I will check whether the very dim (in V-band)
    stars in the V band filter (bolometric magnitude
    is akin to V band magnitude) are low initial mass, low temperature
    stars. The isochrone used to generate those stars will be 10^6.6 years old
    and have solar metallicity"""
    iso1 = synthetic.Isochrone_Binary(6.6, 0, 10,
                                  math.log10(1),
                                  mass_sampling=1)
    
    # Now, consider very dim to be M_V > 17.5.
    # Physical intuition tells us that the low mass stars (e.g.
    # stars closer to M spectral type or brown dwarf area on the
    # HR diagram.
    dim_prims = iso1.primaries[np.where((iso1.primaries['m_ubv_V'] > 17.5))]
    dim_singles = iso1.singles[np.where((iso1.singles['m_ubv_V'] > 17.5))]
    dim_secondaries = iso1.secondaries[np.where((iso1.secondaries['m_ubv_V'] > 17.5))]
    print("Number of primary stars with $M_V < 17.5")
    print(len(dim_prims))
    print("Number of single stars with $M_V < 17.5")
    print(len(dim_singles))
    print("Number of secondary stars with $M_V < 17.5")
    print(len(dim_secondaries))
    assert all(dim_prims['Teff'] < 2000)
    assert all(dim_singles['Teff'] < 2000)
    assert all(dim_secondaries['Teff'] < 2000)
    assert all(dim_prims['mass'] < 1.0)
    assert all(dim_singles['mass'] < 1.0)
    assert all(dim_secondaries['mass'] < 1.0)


def test_cluster_2():
    # Check if the evolution class works fine
    # BPASS Isochrone used to create cluster
    iso1 = synthetic.Isochrone_Binary(9.2, 0.5, 8000,
                                  math.log10(1),
                                  mass_sampling=1)
    # MIST v1 Isochrone used to create cluster
    iso2 = synthetic.IsochronePhot(9.2, 0.5, 8000,
                                   math.log10(1),
                                   recomp=True)
    clus_1 = synthetic.Cluster_w_Binaries(iso1,
                                          imf.imf.IMFSalpeter1955(multiplicity=
                                                          imf.multiplicity.
                                                          MultiplicityUnresolved()),
                                      200000, ifmr=ifmr.IFMR_Spera15())
    clus_2 = synthetic.ResolvedCluster(iso2,
                                       imf.imf.IMFSalpeter1955(multiplicity=imf. multiplicity.
                                                           MultiplicityResolvedDK()),
                                       200000, ifmr=ifmr.IFMR_Spera15())
    assert (np.abs(clus_1.star_systems['systemMass'].sum() -
                  clus_2.star_systems['systemMass'].sum()) ==
            approx(0, rel=10**-2))
    # Now trying this with the custom IMF from Figure 3 of my lunch talk.
    custom_IMF = imf.imf.IMF_broken_powerlaw(np.array([1, 150]),
                                             np.array([-1.7]),
                                             multiplicity=
                                             (imf.multiplicity.
                                              MultiplicityResolvedDK()))
    clus_1 = synthetic.Cluster_w_Binaries(iso1, custom_IMF,
                                      1000000, ifmr=ifmr.IFMR_Spera15())
    clus_2 = synthetic.ResolvedCluster(iso2, custom_IMF,
                                       1000000, ifmr=ifmr.IFMR_Spera15())
    assert np.abs(clus_1['systemMass'].sum() - clus_2['systemMass'].sum()) == approx(0, rel=10**-2)
    # Now trying this with the Kennicutt 1983 from Figure 3 of my lunch talk.
    Kennicutt_IMF = imf.imf.Kennicutt1983(multiplicity = imf.multiplicity.
                                          MultiplicityResolvedDK())
    clus_1 = synthetic.Cluster_w_Binaries(iso1, Kennicutt_IMF,
                                      1000000, ifmr=ifmr.IFMR_Spera15())
    clus_2 = synthetic.ResolvedCluster(iso2, custom_IMF,
                                       1000000, ifmr=ifmr.IFMR_Spera15())
    assert np.abs(clus_1['systemMass'].sum() - clus_2['systemMass'].sum()) == approx(0, rel=10**-2)
    

def test_third_figure_revised():
    """ In this test I will try to make sure that
    BPASS clusters' KLF's do not deviate too widely
    from that of the Observed KLF from the Young Galactic Center Stars
    from Lu et al. 2013. (around 6 million years of age, solar metallicity,
    8000 parsecs from the Earth, and 2.7 AKs extinction).
    For our purposes, we want to make sure that the Observed
    and BPASS counterpart's KLF's are no farther apart than 0.25.
    We may also inspect the KLF of the BPASS and MIST v.1. cluster
    and also the overall initial mass of both BPASS and MIST v1 counterparts
    to the MIST cluster.
    Big Change Made: Use of MultiplicityUnresolved (Lu et. al 2013)
    Use of Merged Models
    IFMR = None is used
    """

    result = load_klf_by_radius(mask_for_log=True)
    magBin = result.Kp[1] - result.Kp[0]
    idx = np.where(result.Kp < 16)[0]
    klf_mag_bins = np.arange(9.0, 17, 1.0)
    binsKp = klf_mag_bins
    binEdges = binsKp[0:-1] + (binsKp[1:] - binsKp[0:-1]) / 2.0


    BPASS_iso = synthetic.Isochrone_Binary(6.78, 2.7,
                                           8000, 0.0,
                                           filters=['nirc2,Kp'])
   # If we want to use control group, let's use multiplicity = None.
    custom_IMF = imf.imf.IMF_broken_powerlaw(np.array([1, 150]),
                                             np.array([-1.7]),
                                             multiplicity=
                                             (imf.multiplicity.
                                              MultiplicityUnresolved()))
    BPASS_Cluster = synthetic.Cluster_w_Binaries(BPASS_iso, custom_IMF,
                                                 1000000,
                                                 ifmr=None)
    
    area = 150.  # arcsec^2
    clus_mass = 1000000
    # Think of how many stars per area would there be in a cluster
    # similar to BPASS_Cluster
    # but with lower mass. We use scaling.
    scaleFactorBPS = ((17000 /
                       clus_mass) /
                      area)
    Merged_iso = synthetic.IsochronePhot(6.78, 2.7, 8000, 0.0,
                                       evo_model=evolution.MergedBaraffePisaEkstromParsec(),
                                       recomp=False,
                                       filters=['nirc2,Kp'])
    Merged_iso = synthetic.IsochronePhot(6.78, 2.7, 8000, 0.0,
                                         recomp=False,
                                         filters=['nirc2,Kp'])
    Mist_Cluster = synthetic.ResolvedCluster(MIST_iso, custom_IMF,
                                             1000000,
                                             ifmr=None)
    Merged_Cluster = synthetic.ResolvedCluster(Merged_iso, custom_IMF,
                                             1000000,
                                             ifmr=None)
    totl_mist =(MIST_Cluster.star_systems['m_nirc2_Kp']
            [np.where((MIST_Cluster.star_systems['isWR'] == 0))[0]])
    totl_merged =(Merged_Cluster.star_systems['m_nirc2_Kp']
                  [np.where((Merged_Cluster.star_systems['isWR'] == 0))[0]])
    totl_BPASS = (BPASS_Cluster.star_systems['m_nirc2_Kp']
               [np.where((~BPASS_Cluster.star_systems['isWR']))[0]])
    mist_scale = scaleFactorBPS
    merged_scale = mist_scale
    weightsMST = np.array([1.0 for x in totl_mist])
    weightsMST *= mist_scale
    weightsBPS = np.array([1.0 for x in totl_BPASS])
    weightsBPS *= scaleFactorBPS
    weightsMerged = np.array([1.0 for x in totl_merged])
    weightsMerged *= merged_scale
    # Binning the K' magnitudes into their respective bins
    # And plotting KLF's of BPASS_Cluster, MIST clusters, and
    # Observed KLF
    # n, bins, patches are exactly the same objects that are
    # described as outputs of matplotlib.hist
    # Right below: BPASS cluster KLF histogram's bin values, bin-edges, and patches
    n, bins, patches = py.hist(totl_mist, bins=binEdges, histtype='step',
                               weights=weightsMST, color='green', label='Merged Model',
                               align='mid', linewidth=1.5)
    # Right below: BPASS cluster KLF histogram's bin values, bin-edges, and patches
    n2, bins2, patches2 = py.hist(totl_BPASS, bins=binEdges, histtype='step',
                                  weights=weightsBPS, color='blue',
                                  label='(BPASS Model KLF)',
                                  align='mid', linewidth=1.5)
    n3, bins3, patches3 = py.hist(totl_merged, bins=binEdges, histtype='step',
                                  weights=weightsMerged, color='blue',
                                  label='(BPASS Model KLF)',
                                  align='mid', linewidth=1.5)
    py.errorbar(result.Kp[idx], result.KLF_ext_cmp_sp_im_noWR[idx],
                fmt='ro-', xerr=magBin/2.0, capsize=0, linewidth=2)
    py.errorbar(result.Kp[idx] ,result.KLF_ext_cmp_sp_im_noWR[idx],
                fmt='ro-', yerr=result.eKLF_ext_cmp_sp_im_noWR[idx],
                linewidth=2,
                label='Observed')
    py.xlim(8.5, 15.5)
    py.xlabel('Kp magnitude')
    py.ylabel("stars / (arcsecond^2 mag)")
    py.title("KLF's at Age = %d Myr" % (10**(6.78 - 6)), fontsize=14)
    py.legend(loc='upper left', numpoints=1)
    py.savefig('Comparisons_w_Real_Data.png')
    # Now compare differences in the heights of BPASS cluster KLF and given KLF
    for x in range(len(binEdges)):
        assert np.abs(result.KLF_ext_cmp_sp_im_noWR[idx][x] - n[x]) == approx(0, abs=0.25)
    print("BPASS IMF does not seem too far off when compared to the BPASS IMF")
    # Now compare differences in the BPASS and MIST KLF's (or values
    # of bins in histograms for both KLF's)
    for x in range(len(binEdges)):
        assert np.abs(n2[x] - n[x]) == approx(0, abs=0.125)
        assert np.abs(n2[x] - n3[x]) == approx(0, abs=0.125)
    # Now we try to generate a figure closer to the Figure 1 of Lu et al 2013
    custom_IMF2 = 
    BPASS_Cluster = synthetic.Cluster_w_Binaries(BPASS_iso, custom_IMF,
                                                 1000000,
                                                 ifmr=None)
    
def test_third_figure_w_Tab2Sol1():
    """ In this test I will try to make sure that
    BPASS clusters' KLF's do not deviate too widely
    from that of the Observed KLF from the Young Galactic Center Stars
    from Lu et al. 2013. (around 6 million years of age, solar metallicity,
    8000 parsecs from the Earth, and 2.7 AKs extinction).
    For our purposes, we want to make sure that the Observed
    and BPASS counterpart's KLF's are no farther apart than 0.25.
    We may also inspect the KLF of the BPASS and MIST v.1. cluster
    and also the overall initial mass of both BPASS and MIST v1 counterparts
    to the MIST cluster.
    Big Change Made: Use of MultiplicityUnresolved (Lu et. al 2013)
    Use of Merged Models
    IFMR = None is used
    Using Table 2 Solution 1 from Lu et. al 2013
    """

    result = load_klf_by_radius(mask_for_log=True)
    magBin = result.Kp[1] - result.Kp[0]
    idx = np.where(result.Kp < 16)[0]
    klf_mag_bins = np.arange(9.0, 17, 1.0)
    binsKp = klf_mag_bins
    binEdges = binsKp[0:-1] + (binsKp[1:] - binsKp[0:-1]) / 2.0


    BPASS_iso = synthetic.Isochrone_Binary(6.78, 2.7,
                                           8000, 0.0,
                                           filters=['nirc2,Kp'])
   # If we want to use control group, let's use multiplicity = None.
    custom_IMF = imf.imf.IMF_broken_powerlaw(np.array([1, 150]),
                                             np.array([-1.7]),
                                             multiplicity=
                                             (imf.multiplicity.
                                              MultiplicityUnresolved()))
    BPASS_Cluster = synthetic.Cluster_w_Binaries(BPASS_iso, custom_IMF,
                                                 1000000,
                                                 ifmr=None)
    
    area = 150.  # arcsec^2
    clus_mass = 1000000
    # Think of how many stars per area would there be in a cluster
    # similar to BPASS_Cluster
    # but with lower mass. We use scaling.
    scaleFactorBPS = ((17000 /
                       clus_mass) /
                      area)
    Merged_iso = synthetic.IsochronePhot(6.62, 2.7, 8000, 0.0,
                                       evo_model=evolution.MergedBaraffePisaEkstromParsec(),
                                       recomp=False,
                                       filters=['nirc2,Kp'])
    Merged_iso = synthetic.IsochronePhot(6.62, 2.7, 8000, 0.0,
                                         recomp=False,
                                         filters=['nirc2,Kp'])
    Mist_Cluster = synthetic.ResolvedCluster(MIST_iso, custom_IMF,
                                             1000000,
                                             ifmr=None)
    Merged_Cluster = synthetic.ResolvedCluster(Merged_iso, custom_IMF,
                                             1000000,
                                             ifmr=None)
    totl_mist =(MIST_Cluster.star_systems['m_nirc2_Kp']
            [np.where((MIST_Cluster.star_systems['isWR'] == 0))[0]])
    totl_merged =(Merged_Cluster.star_systems['m_nirc2_Kp']
                  [np.where((Merged_Cluster.star_systems['isWR'] == 0))[0]])
    totl_BPASS = (BPASS_Cluster.star_systems['m_nirc2_Kp']
               [np.where((~BPASS_Cluster.star_systems['isWR']))[0]])
    mist_scale = scaleFactorBPS
    merged_scale = mist_scale
    weightsMST = np.array([1.0 for x in totl_mist])
    weightsMST *= mist_scale
    weightsBPS = np.array([1.0 for x in totl_BPASS])
    weightsBPS *= scaleFactorBPS
    weightsMerged = np.array([1.0 for x in totl_merged])
    weightsMerged *= merged_scale
    # Binning the K' magnitudes into their respective bins
    # And plotting KLF's of BPASS_Cluster, MIST clusters, and
    # Observed KLF
    # n, bins, patches are exactly the same objects that are
    # described as outputs of matplotlib.hist
    # Right below: BPASS cluster KLF histogram's bin values, bin-edges, and patches
    n, bins, patches = py.hist(totl_mist, bins=binEdges, histtype='step',
                               weights=weightsMST, color='green', label='Merged Model',
                               align='mid', linewidth=1.5)
    # Right below: BPASS cluster KLF histogram's bin values, bin-edges, and patches
    n2, bins2, patches2 = py.hist(totl_BPASS, bins=binEdges, histtype='step',
                                  weights=weightsBPS, color='blue',
                                  label='(BPASS Model KLF)',
                                  align='mid', linewidth=1.5)
    n3, bins3, patches3 = py.hist(totl_merged, bins=binEdges, histtype='step',
                                  weights=weightsMerged, color='blue',
                                  label='(BPASS Model KLF)',
                                  align='mid', linewidth=1.5)
    py.errorbar(result.Kp[idx], result.KLF_ext_cmp_sp_im_noWR[idx],
                fmt='ro-', xerr=magBin/2.0, capsize=0, linewidth=2)
    py.errorbar(result.Kp[idx] ,result.KLF_ext_cmp_sp_im_noWR[idx],
                fmt='ro-', yerr=result.eKLF_ext_cmp_sp_im_noWR[idx],
                linewidth=2,
                label='Observed')
    py.xlim(8.5, 15.5)
    py.xlabel('Kp magnitude')
    py.ylabel("stars / (arcsecond^2 mag)")
    py.title("KLF's at Age = %d Myr" % (10**(6.62 - 6)), fontsize=14)
    py.legend(loc='upper left', numpoints=1)
    py.savefig('Comparisons_w_Real_Data.png')
    # Now compare differences in the heights of BPASS cluster KLF and given KLF
    for x in range(len(binEdges)):
        assert np.abs(result.KLF_ext_cmp_sp_im_noWR[idx][x] - n[x]) == approx(0, abs=0.25)
    print("BPASS IMF does not seem too far off when compared to the BPASS IMF")
    # Now compare differences in the BPASS and MIST KLF's (or values
    # of bins in histograms for both KLF's)
    for x in range(len(binEdges)):
        assert np.abs(n2[x] - n[x]) == approx(0, abs=0.1)
        assert np.abs(n2[x] - n3[x]) == approx(0, abs=0.1)
    
