##########################################################
#
#This IFMF comes from Raithel et al. 2017
#https://arxiv.org/pdf/1712.00021.pdf
#
#########################################################

import numpy as np

class IFMF(object):
    def __init__(self):
        """
        The IFMF class.
        """

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
        return f_ej * BH_mass_core_low(MZAMS) + (1 - f_ej) * BH_mass_all(MZAMS)

    def NS_mass(self, MZAMS):
        """
        Paper: 9 < MZAMS 120
        Simplify to just return one value
        """
        return 1.6

    def WD_mass(self, MZAMS):
        """
        From Kalirai+07
        1.16 < MZAMS < 6.5
        FIXME: need to extend these ranges...
        """
        return 0.109*MZAMS + 0.394
    
    def generate_death_mass_distribution(self, mass_array):
        """
        WD: typecode = 1
        NS: typecode = 2
        BH: typecode = 3
        MZAMS > 120: big
        MZAMS < 0.5: small
        """
        for i in np.arange(len(mass_array)):
            MZAMS = mass_array[i]
            n = random.randint(1,101)
            if (MZAMS >= 0.5) and (MZAMS < 9):
                typecode = 1
                return WD_mass(MZAMS), typecode
            elif (MZAMS >=9) and (MZAMS <= 15):
                typecode = 2
                return NS_mass(MZAMS), typecode
            elif (MZAMS > 15) and (MZAMS <= 17.8):
                if n > 68:
                    typecode = 3
                    return BH_mass_low(MZAMS, 0.9), typecode
                else:
                    typecode = 2
                    return NS_mass(MZAMS), typecode
            elif (MZAMS > 17.8) and (MZAMS <= 18.5):
                if n > 83:
                    typecode = 3
                    return BH_mass_low(MZAMS, 0.9), typecode
                else:
                    typecode = 2
                    return NS_mass(MZAMS), typecode
            elif (MZAMS > 18.5) and (MZAMS <= 21.7):
                if n > 50:
                    typecode = 3
                    return BH_mass_low(MZAMS, 0.9), typecode
                else:
                    typecode = 2
                    return NS_mass(MZAMS), typecode
            elif (MZAMS > 21.7) and (MZAMS <= 25.2):
                typecode = 3
                return BH_mass_low(MZAMS, 0.9), typecode
            elif (MZAMS > 25.2) and (MZAMS <= 27.5):
                if n > 65:
                    typecode = 3
                    return BH_mass_low(MZAMS, 0.9), typecode
                else: 
                    typecode = 2
                    return NS_mass(MZAMS), typecode
            elif (MZAMS > 27.5) and (MZAMS <= 60):
                typecode = 3
                if MZAMS > 42.22:
                    return BH_mass_high(MZAMS), typecode
                else:
                    return BH_mass_low(MZAMS, 0.9), typecode
            elif (MZAMS > 60) and (MZAMS <= 120):
                if n > 40:
                    typecode = 3
                    return BH_mass_high(MZAMS), typecode
                else:
                    typecode = 2
                    return NS_mass(MZAMS), typecode
            elif MZAMS < 0.5:
                typecode = 'small'
                return 0, typecode
            else:
                typecode = 'big'
                return 0, typecode
                    
    
                                                                                       
