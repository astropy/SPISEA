# Test functions for the different stellar evolution and atmosphere models
import pdb

def test_evolution_models():
    """
    Test to make sure the different evolution models work
    """
    from popstar import evolution

    # Age ranges to test
    age_young_arr = [6.7, 7.9]
    age_all_arr = [6.7, 8.0, 9.75]

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
        print('Done {0}'.format(evo))
        
    return

def test_atmosphere_models():
    """
    Test the rebinned atmosphere models used for synthetic photometry
    """
    from popstar import atmospheres as atm

    # Array of atmospheres
    atm_arr = [atm.get_merged_atmosphere, atm.get_castelli_atmosphere, atm.get_phoenixv16_atmosphere, atm.get_BTSettl_2015_atmosphere,
                   atm.get_kurucz_atmosphere, atm.get_phoenix_atmosphere]

    # Array of metallicities
    metals_range = [-2.0, 0, 0.15]
    metals_solar = [0]
    metals_arr = [metals_solar, metals_range, metals_range, metals_solar, metals_range, metals_range]

    assert len(atm_arr) == len(metals_arr)

    # Loop through models, testing if them work
    for ii in range(len(atm_arr)):
        atm = atm_arr[ii]

        # Loop through metallicities
        for jj in metals_arr[ii]:
            try:
                test = atm(metallicity=jj)
            except:
                print('TEST FAILED: {0}, metal = {1}'.format(atm, jj))
                
        print('Done {0}'.format(atm))

    return
