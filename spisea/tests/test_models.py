# Test functions for the different stellar evolution and atmosphere models
import numpy as np
import pdb

def test_evo_model_grid_num():
    """
    Make sure evolution models have both evo_grid_num 
    and evo_grid_min (e.g., make sure these functions
    are working). Try it on one evolution model here;
    we'll test on all evo models in another function.
    """
    from spisea import evolution
    
    # Make MIST evolution model, check evo grid variables
    evo = evolution.MISTv1()
    assert isinstance(evo.evo_grid_min, float)
    
    return

def test_evolution_models():
    """
    Test to make sure the different evolution models work
    """
    from spisea import evolution

    # Age ranges to test
    age_young_arr = [6.7, 7.9]
    age_all_arr = [6.7, 8.0, 9.7]
    age_all_MIST_arr = [5.2, 6.7, 9.7, 10.13]
    bd_young_test = [6.0, 6.5, 7.4]

    # Metallicity ranges to test (if applicable)
    metal_range = [-2.5, -1.5, 0, 0.25, 0.4]
    metal_solar = [0]
    metal_Marley = [-0.5, 0.0, 0.5]

    # Array of evolution models to test
    evo_models = [evolution.MISTv1(version=1.2), evolution.MergedBaraffePisaEkstromParsec(), 
                      evolution.Parsec(), evolution.Baraffe15(), evolution.Ekstrom12(), evolution.Pisa(), 
                      evolution.Phillips2020(), evolution.Marley2021(), 
                      evolution.MergedPhillipsBaraffePisaEkstromParsec()]

    
    # Array of age_ranges for the specific evolution models to test
    age_vals = [age_all_MIST_arr, age_all_arr, age_all_arr, age_young_arr, age_young_arr, age_young_arr, age_all_arr, age_all_arr, bd_young_test]

    # Array of metallicities for the specific evolution models to test
    metal_vals = [metal_range, metal_solar, metal_solar, metal_solar, metal_solar, metal_solar, metal_solar, metal_Marley, metal_solar]

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

                    # Make sure the actual isochrone metallicity taken is
                    # indeed the closest to the desired metallicity (e.g., closest point
                    # in evo grid)
                    z_ratio = np.log10(np.array(evo.z_list) / evo.z_solar)
                    closest_idx = np.where( abs(z_ratio - jj) == min(abs(z_ratio - jj)) )[0][0]
                    expected_val = z_ratio[closest_idx]

                    assert np.isclose(test.meta['metallicity_act'], expected_val, atol=0.01)

                except:
                    raise Exception('EVO TEST FAILED: {0}, age = {1}, metal = {2}'.format(evo, kk, jj))

        print('Done {0}'.format(evo))
        
    return

def test_atmosphere_models():
    """
    Test the rebinned atmosphere models used for synthetic photometry
    """
    from spisea import atmospheres as atm

    # Array of atmospheres
    atm_arr = [atm.get_merged_atmosphere, atm.get_castelli_atmosphere, atm.get_phoenixv16_atmosphere, 
               atm.get_BTSettl_2015_atmosphere, atm.get_BTSettl_atmosphere, atm.get_kurucz_atmosphere, 
               atm.get_phoenix_atmosphere, atm.get_wdKoester_atmosphere, atm.get_Phillips2020_atmosphere, 
               atm.get_Meisner2023_atmosphere]

    # Array of metallicities
    metals_range = [-2.0, 0, 0.15]
    bd_metals_range = [-1.0, -0.5, 0, 0.3]
    metals_solar = [0]
    metals_arr = [metals_solar, metals_range, metals_range, metals_solar, metals_range, metals_range, metals_range, 
                  metals_solar, metals_solar, bd_metals_range]

    assert len(atm_arr) == len(metals_arr)

    # Loop through models, testing if them work
    for ii in range(len(atm_arr)):
        atm_func = atm_arr[ii]

        # Loop through metallicities
        for jj in metals_arr[ii]:
            try:
                test = atm_func(metallicity=jj)
            except:
                raise Exception('ATM TEST FAILED: {0}, metal = {1}'.format(atm_func, jj))
                
        print('Done {0}'.format(atm_func))
        
    # Test get_merged_atmospheres at different temps
    temp_range = [200, 1000, 2000, 3500, 4000, 5250, 6000, 12000]
    atm_func = atm.get_merged_atmosphere
    for ii in bd_metals_range:
        for jj in temp_range:
            try:
                test = atm_func(metallicity=ii, temperature=jj, verbose=True)
            except:
                raise Exception('ATM TEST FAILED: {0}, metal = {1}, temp = {2}'.format(atm_func, ii, jj))


    print('get_merged_atmosphere: all temps/metallicities passed')
    
    # Test get_bb_atmosphere at different temps
    # This func only requests temp
    temp_range = [1000, 2000, 3500, 4000, 5250, 6000, 12000]
    atm_func = atm.get_bb_atmosphere
    for jj in temp_range:
        try:
            test = atm_func(temperature=jj, verbose=True)
        except:
            raise Exception('ATM TEST FAILED: {0}, temp = {2}'.format(atm_func, jj))
    
    print('get_bb_atmosphere: all temps passed')

    # Test get_bd_atmosphere at different temps
    # This func only requests temp
    temp_range = [250, 400, 500, 750, 950, 1200]
    atm_func = atm.get_bd_atmosphere
    for jj in temp_range:
        try:
            test = atm_func(temperature=jj, verbose=True)
        except:
            raise Exception('ATM TEST FAILED: {0}, temp = {1}'.format(atm_func, jj))
    
    print('get_bd_atmosphere: all temps passed')
    
    return

def test_filters():
    """
    Test to make sure all of the filters work as expected
    """
    from spisea import synthetic

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
                     'naco,IB_2.00','naco,IB_2.03', 'naco,IB_2.06', 'naco,IB_2.24',
                     'naco,IB_2.27', 'naco,IB_2.30', 'naco,IB_2.33', 'naco,IB_2.36',
                     'nirc1,K', 'nirc1,H', 'nirc2,J', 'nirc2,H',
                     'nirc2,Kp', 'nirc2,K', 'nirc2,Lp', 'nirc2,Hcont',
                     'nirc2,FeII', 'nirc2,Brgamma', 'ps1,z',
                     'ps1,g', 'ps1,r','ps1,i', 'ps1,y',
                     'ukirt,J', 'ukirt,H', 'ukirt,K',
                     'vista,Y', 'vista,Z', 'vista,J',
                     'vista,H',  'vista,Ks', 'ztf,g', 'ztf,r', 'ztf,i',
                     'hawki,J', 'hawki,H', 'hawki,Ks', 'roman,wfi,f062',
                     'roman,wfi,f087', 'roman,wfi,f106', 'roman,wfi,f129',
                     'roman,wfi,f158', 'roman,wfi,w146', 'roman,wfi,f213',
                     'roman,wfi,f184', 'rubin,g', 'rubin,i', 'rubin,r',
                     'rubin,u', 'rubin,z', 'rubin,y']

    # Loop through filters to test that they work: get_filter_info
    for ii in filt_list:
        try:
            filt = synthetic.get_filter_info(ii, rebin=True, vega=vega)
        except:
            raise Exception('get_filter_info TEST FAILED for {0}'.format(ii))

    print('get_filter_info pass')
    
    # Loop through filters to test that they work: get_obs_str
    for ii in filt_list:
        try:
            # Test going from col_name to obs_str
            col_name = synthetic.get_filter_col_name(ii)
            obs_str = synthetic.get_obs_str('m_{0}'.format(col_name))
            # Does the obs_str work?
            filt_info = synthetic.get_filter_info(obs_str)
        except:
            raise Exception('get_obs_str TEST FAILED for {0}'.format(ii)) 
            
    print('get_obs_str pass')
    print('Filters done')

    return
