##########################################################
#
#
#This IFMR comes from Raithel et al. 2017
#https://arxiv.org/pdf/1712.00021.pdf
#
#
#########################################################

import numpy as np


class IFMR(object):
    """
    The IFMR base class. The IFMR is a combination of the 
    WD IFMR from 
    `Kalirai et al. (2008) <https://ui.adsabs.harvard.edu/abs/2008ApJ...676..594K/abstract>`_
    and the NS/BH IFMR from
    `Raithel et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018ApJ...856...35R/abstract>`_.

    See Lam et al. (submitted) for more details.
    """
    def __init__(self):
        pass

    def BH_mass_core_low(self, MZAMS):
        """                                                                                                      
        Eqn (1)                                                                                                  
        Paper: 15 < MZAMS < 40                                                                                   
        Us extending: 15 < MZAMS < 42.22                                                                         
        """
        return -2.024 + 0.4130*MZAMS

    def BH_mass_all_low(self, MZAMS):
        """                                                                                                      
        Eqn (2)                                                                                                  
        Paper: 15 < MZAMS < 40                                                                                   
        Us extending: 15 < MZAMS < 42.22                                                                         
        """
        return 16.28 + 0.00694 * (MZAMS - 21.872) - 0.05973 * (MZAMS - 21.872)**2 + 0.003112 * (MZAMS - 21.872)**3

    def BH_mass_high(self, MZAMS):
        """                                                                                                      
        Eqn (3)                                                                                                  
        Paper: 45 < MZAMS < 120                                                                                  
        Us extending: 42.22 < MZAMS < 120                                                                        
        """
        return 5.795 + 1.007 * 10**9 * MZAMS**-4.926

    def BH_mass_low(self, MZAMS, f_ej):
        """                                                                                                      
        Eqn (4)                                                                                                  
        Paper: 15 < MZAMS < 40                                                                                   
        Us extending: 15 < MZAMS < 42.22                                                                         
        """
        return f_ej * self.BH_mass_core_low(MZAMS) + (1 - f_ej) * self.BH_mass_all_low(MZAMS)

    def NS_mass(self, MZAMS):
        """                                                                                                      
        Paper: 9 < MZAMS 120                                                                                     
        Simplify to just return one value                                                                        
        """
        return 1.6 * np.ones(len(MZAMS))

    def WD_mass(self, MZAMS):
        """                                                                                                      
        From Kalirai+07                                                                                          
        1.16 < MZAMS < 6.5                                                                                       
        FIXME: need to extend these ranges...                                                                    
        """
        return 0.109*MZAMS + 0.394

    def generate_death_mass(self, mass_array):
        """
        The top-level function that assigns the remnant type 
        and mass based on the stellar initial mass. 
        
        Parameters
        ----------
        mass_arr: array of floats
            Array of initial stellar masses. Units are
            M_sun.

        Notes
        ------
        The output typecode tells what compact object formed:
        
        * WD: typecode = 101
        * NS: typecode = 102
        * BH: typecode = 103

        A typecode of value -1 means you're outside the range of 
        validity for applying the ifmr formula. 

        A remnant mass of -99 means you're outside the range of 
        validity for applying the ifmr formula.

        Range of validity: 0.5 < MZAMS < 120

        Returns
        -------
        output_arr: 2-element array
            output_array[0] contains the remnant mass, and 
            output_array[1] contains the typecode
        """

        #output_array[0] holds the remnant mass
        #output_array[1] holds the remnant type
        output_array = np.zeros((2, len(mass_array)))

        #Random array to get probabilities for what type of object will form
        random_array = np.random.randint(1, 1001, size = len(mass_array))

        codes = {'WD': 101, 'NS': 102, 'BH': 103}
        
        """
        The id_arrays are to separate all the different formation regimes
        """
        id_array0 = np.where((mass_array < 0.5) | (mass_array >= 120))
        output_array[0][id_array0] = -99 * np.ones(len(id_array0))
        output_array[1][id_array0]  = -1 * np.ones(len(id_array0))

        id_array1 = np.where((mass_array >= 0.5) & (mass_array < 9))
        output_array[0][id_array1] = self.WD_mass(mass_array[id_array1])
        output_array[1][id_array1]= codes['WD']

        id_array2 = np.where((mass_array >= 9) & (mass_array < 15))
        output_array[0][id_array2] = self.NS_mass(mass_array[id_array2])
        output_array[1][id_array2] = codes['NS']

        id_array3_BH = np.where((mass_array >= 15) & (mass_array < 17.8) & (random_array > 679))
        output_array[0][id_array3_BH] = self.BH_mass_low(mass_array[id_array3_BH], 0.9)
        output_array[1][id_array3_BH] = codes['BH']

        id_array3_NS = np.where((mass_array >= 15) & (mass_array < 17.8) & (random_array <= 679))
        output_array[0][id_array3_NS] = self.NS_mass(mass_array[id_array3_NS])
        output_array[1][id_array3_NS] = codes['NS']

        id_array4_BH = np.where((mass_array >= 17.8) & (mass_array < 18.5) & (random_array > 833))
        output_array[0][id_array4_BH]= self.BH_mass_low(mass_array[id_array4_BH], 0.9)
        output_array[1][id_array4_BH] = codes['BH']
        
        id_array4_NS = np.where((mass_array >= 17.8) & (mass_array < 18.5) & (random_array <= 833))
        output_array[0][id_array4_NS] = self.NS_mass(mass_array[id_array4_NS])
        output_array[1][id_array4_NS] = codes['NS']

        id_array5_BH = np.where((mass_array >= 18.5) & (mass_array < 21.7) & (random_array > 500))
        output_array[0][id_array5_BH] = self.BH_mass_low(mass_array[id_array5_BH], 0.9)
        output_array[1][id_array5_BH] = codes['BH']
        
        id_array5_NS = np.where((mass_array >= 18.5) & (mass_array < 21.7) & (random_array <= 500))
        output_array[0][id_array5_NS] = self.NS_mass(mass_array[id_array5_NS])
        output_array[1][id_array5_NS] = codes['NS']

        id_array6 = np.where((mass_array >= 21.7) & (mass_array < 25.2))
        output_array[0][id_array6] = self.BH_mass_low(mass_array[id_array6], 0.9)
        output_array[1][id_array6]= codes['BH']

        id_array7_BH = np.where((mass_array >= 25.2) & (mass_array < 27.5) & (random_array > 652))
        output_array[0][id_array7_BH] = self.BH_mass_low(mass_array[id_array7_BH], 0.9)
        output_array[1][id_array7_BH] = codes['BH']
        
        id_array7_NS = np.where((mass_array >= 25.2) & (mass_array < 27.5) & (random_array <= 652))
        output_array[0][id_array7_NS] = self.NS_mass(mass_array[id_array7_NS])
        output_array[1][id_array7_NS] = codes['NS']

        id_array8 = np.where((mass_array >= 27.5) & (mass_array < 42.22))
        output_array[0][id_array8] = self.BH_mass_low(mass_array[id_array8], 0.9)
        output_array[1][id_array8] = codes['BH']

        id_array9 = np.where((mass_array >= 42.22) & (mass_array < 60))
        output_array[0][id_array9] = self.BH_mass_high(mass_array[id_array9])
        output_array[1][id_array9] = codes['BH']

        id_array10_BH = np.where((mass_array >= 60) & (mass_array < 120) & (random_array > 400))
        output_array[0][id_array10_BH] = self.BH_mass_high(mass_array[id_array10_BH])
        output_array[1][id_array10_BH] = codes['BH']
        
        id_array10_NS = np.where((mass_array >= 60) & (mass_array < 120) & (random_array <= 400))
        output_array[0][id_array10_NS] = self.NS_mass(mass_array[id_array10_NS])
        output_array[1][id_array10_NS] = codes['NS']

        return(output_array)
