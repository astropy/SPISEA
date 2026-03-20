from spisea import evolution, exceptions
import copy
import os
import pdb


def test_grid_number_exception():
    """
    Make sure check_evo_grid_number works as expected
    """
    #----Case 1: installed model grid is too low (raise error)----#
    # Get installed grid num
    models_dir = '{0}/evolution/'.format(os.environ['SPISEA_MODELS'])
    installed_grid = evolution.get_installed_grid_num(models_dir)

    # Make sure required grid is too large
    required_grid = installed_grid + 1

    # Make sure correct error is raised, otherwise fail
    try:
        evolution.check_evo_grid_number(required_grid, models_dir)
    except exceptions.ModelMismatch:
        pass
    else:
        raise Exception('check_evo_grid_number failed')

    # Case 2: installed model grid is same (no error)
    required_grid = copy.deepcopy(installed_grid)

    evolution.check_evo_grid_number(required_grid, models_dir)

    # Case 3: installed model grid is higher than required grid (no error)
    required_grid = installed_grid - 1.0

    evolution.check_evo_grid_number(required_grid, models_dir)    

    return