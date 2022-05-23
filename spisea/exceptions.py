import pdb


class ModelMismatch(Exception):
    """
    Error raised when installed model grid version number is mismatched
    with the user's version of SPISEA code
    """
    def __init__(self, required_num, grid_num, model_type):
        assert (model_type == 'evolution') | (model_type == 'atmosphere')
        
        if model_type == 'evolution':
            model_file = 'spisea_models.tar.gz'
        elif model_type == 'atmosphere':
            model_file = 'spisea_cdbs.tar.gz'
        
        # Compose error message
        str1 = 'WARNING: Desired {0} model requires model grid version >= {1},'.format(model_type, required_num)
        str2 = 'but model grid version {0} is installed.'.format(grid_num)
        str3 = 'Please re-download {0} from installation documentation to resolve mismatch.'.format(model_file)
        self.message = '{0} {1} {2}'.format(str1, str2, str3)

        # Define error
        super().__init__(self.message)
