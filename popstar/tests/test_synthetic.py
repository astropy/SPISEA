def test_isochrone():
    from popstar import synthetic as syn

    logAge = 6.7
    AKs = 2.7
    distance = 4000
    
    iso = syn.Isochrone(logAge, AKs, distance)

    assert iso.points.meta['LogAge'] == logAge
    assert iso.points.meta['AKs'] == AKs
    assert iso.points.meta['distance'] == distance
    assert len(iso.points) > 100

    return

def test_cluster():

    # multi = multiplicity.MultiplicityUnresolved()
    # imf = imf.Kroupa_2001(multiplicty=multi)
    # evo = evolution.MergedPisaEkstromParsec()
    # atm_func = atm.get_merged_atmosphere

    # log_age = 6.5
    # AKs = 1.0
    # distance = 8000.0

    # if resolved:
    #     cluster = ResolvedCluster(log_age, AKs, distance, imf, evo, atm_func)
    # else:
    #     cluster = UnresolvedCluster(log_age, AKs, distance, imf, evo, atm_func)

    # # Plot the spectrum of the most massive star
    # idx = cluster.mass.argmax()
    # plt.clf()
    # plt.plot(cluster.stars[idx].wave, cluster.stars[idx].flux, 'k.')

    # # Plot an integrated spectrum of the whole cluster.
    # wave, flux = cluster.get_integrated_spectrum()
    # plt.clf()
    # plt.plot(wave, flux, 'k.')

    return
    
