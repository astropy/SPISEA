# Test functions for the different stellar evolution and atmosphere models
import pdb

def test_evolution_models():
    """
    Test to make sure the different evolution models work
    """
    from popstar import evolution

    # Age ranges to test
    age_young_arr = [6.7, 7.9]
    age_all_arr = [6.7, 8.0, 9.7]

    # Metallicity ranges to test (if applicable)
    metal_range = [-2.5, 0, 0.4]
    metal_solar = [0]

    # Array of evolution models to test
    evo_models = [evolution.MISTv1(version=1.2), evolution.MergedBaraffePisaEkstromParsec(), evolution.MergedSiessGenevaPadova(),
                      evolution.Parsec(), evolution.Baraffe15(), evolution.Ekstrom12(), evolution.Pisa()]

    
    # Array of age_ranges for the specific evolution models to test
    age_vals = [age_all_arr, age_all_arr, age_all_arr, age_all_arr, age_young_arr, age_young_arr, age_young_arr]

    # Array of metallicities for the specific evolution models to test
    metal_vals = [metal_range, metal_solar, metal_solar, metal_solar, metal_solar, metal_solar, metal_solar]

    assert len(evo_models) == len(age_vals) == len(metal_vals)

    # Loop through models, testing if them work
    for ii in range(len(evo_models)):
        evo = evo_models[ii]

        # Loop through metallicities
        for jj in metal_vals[ii]:
            # Loop through ages
            for kk in age_vals[ii]:
                try:
                    test = evo.isochrone(age=10**kk, metallicity=jj)
                except:
                    print('TEST FAILED: {0}, age = {1}, metal = {2}'.format(evo, kk, jj))
                    pdb.set_trace()
        print('Done {0}'.format(evo))
        
    return

def test_atmosphere_models():
    """
    Test the rebinned atmosphere models used for synthetic photometry
    """
    from popstar import atmospheres as atm

    # Array of atmospheres
    atm_arr = [atm.get_merged_atmosphere, atm.get_castelli_atmosphere, atm.get_phoenixv16_atmosphere, atm.get_BTSettl_2015_atmosphere,
                   atm.get_BTSettl_atmosphere, atm.get_kurucz_atmosphere, atm.get_phoenix_atmosphere]

    # Array of metallicities
    metals_range = [-2.0, 0, 0.15]
    metals_solar = [0]
    metals_arr = [metals_solar, metals_range, metals_range, metals_solar, metals_range, metals_range, metals_range]

    assert len(atm_arr) == len(metals_arr)

    # Loop through models, testing if them work
    for ii in range(len(atm_arr)):
        atm_func = atm_arr[ii]

        # Loop through metallicities
        for jj in metals_arr[ii]:
            try:
                test = atm_func(metallicity=jj)
            except:
                print('TEST FAILED: {0}, metal = {1}'.format(atm_func, jj))
                
        print('Done {0}'.format(atm_func))
        
    # One last test: get_merged_atmospheres at different temps
    temp_range = [2000, 4000, 6000, 12000]
    atm = atm.get_merged_atmosphere
    for ii in metal_range:
        for jj in temp_range:
            try:
                test = atm_func(metallicity=ii, temperature=jj)
            except:
                print('TEST FAILED: {0}, metal = {1}, temp = {2}'.format(atm_func, ii))            

    return

def test_filters():
    """
    Test to make sure all of the filters work as expected
    """
    from popstar import synthetic

    # Define vega spectrum
    vega = synthetic.Vega()
    
    # Filter list to test
    filt_list = ['wfc3,ir,f127m','acs,wfc1,f814w',
                     '2mass,J', '2mass,H','2mass,Ks',
                     'ctio_osiris,K', 'ctio_osiris,H',
                     'ubv,U', 'ubv,B', 'ubv,V', 'ubv,R',
                     'ubv,I', 'jg,J', 'jg,H', 'jg,K',
                     'decam,y', 'decam,i', 'decam,z',
                     'decam,u', 'decam,g', 'decam,r',
                     'gaia,dr2_rev,G', 'gaia,dr2_rev,Gbp', 'gaia,dr2_rev,Grp',
                     'jwst,F070W', 'jwst,F090W', 'jwst,F115W', 'jwst,F140M',
                     'jwst,F150W', 'jwst,F150W2', 'jwst,F162M', 'jwst,F164N',
                     'jwst,F182M', 'jwst,F187N', 'jwst,F200W', 'jwst,F212N',
                     'jwst,F210M','jwst,F250M', 'jwst,F277W', 'jwst,F300M',
                     'jwst,F322W2', 'jwst,F323N', 'jwst,F335M', 'jwst,F356W',
                     'jwst,F360M', 'jwst,F405N', 'jwst,F410M', 'jwst,F430M',
                     'jwst,F444W', 'jwst,F460M', 'jwst,F466N', 'jwst,F470N',
                     'jwst,F480M', 'naco,J', 'naco,H', 'naco,Ks',
                     'nirc1,K', 'nirc1,H', 'nirc2,J', 'nirc2,H',
                     'nirc2,Kp', 'nirc2,K', 'nirc2,Lp', 'nirc2,Hcont',
                     'nirc2,FeII', 'nirc2,Brgamma', 'ps1,z',
                     'ps1,g', 'ps1,r','ps1,i', 'ps1,y',
                     'ukirt,J', 'ukirt,H', 'ukirt,K',
                     'vista,Y', 'vista,Z', 'vista,J',
                     'vista,H',  'vista,Ks']

    # Loop through filters to test that they work
    for ii in filt_list:
        try:
            filt = synthetic.get_filter_info(ii, rebin=True, vega=vega)
        except:
            print('TEST FAILED for {0}'.format(ii))

    print('Filters done')

    return
