##########################################################
#
#
#This IFMF comes from Raithel et al. 2017
#https://arxiv.org/pdf/1712.00021.pdf
#
#
#Casey's note to self on how to run this:
#from ifmf import IFMF
#t = IFMF()
#t.generate_death_mass_distribution([array here])
#
#
#########################################################

import numpy as np

class IFMF(object):
    def __init__(self):
        pass
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
        return f_ej * self.BH_mass_core_low(MZAMS) + (1 - f_ej) * self.BH_mass_all_low(MZAMS)

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
        typecode tells what compact object formed:
        WD: typecode = 1
        NS: typecode = 2
        BH: typecode = 3
        typecode of value -1 means you're outside the range of validity for applying the ifmf formula
        remnant mass of -99 means you're outside the range of validity for applying the ifmf formula
        range of validity: 0.5 < MZAMS < 120
        output_array[0] contains the remnant mass
        output_array[1] contains the typecode
        
        """
        output_array = np.zeros((2, len(mass_array)))
        for i in np.arange(len(mass_array)):
            MZAMS = mass_array[i]
            n = np.random.randint(1,101)
            if (MZAMS >= 0.5) and (MZAMS < 9):
                output_array[0][i] = self.WD_mass(MZAMS)
                output_array[1][i] = 1
            elif (MZAMS >=9) and (MZAMS <= 15):
                output_array[0][i] = self.NS_mass(MZAMS)
                output_array[1][i] = 2
            elif (MZAMS > 15) and (MZAMS <= 17.8):
                if n > 68:
                    output_array[0][i] = self.BH_mass_low(MZAMS, 0.9)
                    output_array[1][i] = 3
                else:
                    output_array[0][i] = self.NS_mass(MZAMS)
                    output_array[1][i] = 2
            elif (MZAMS > 17.8) and (MZAMS <= 18.5):
                if n > 83:
                    output_array[0][i] = self.BH_mass_low(MZAMS, 0.9)
                    output_array[1][i] = 3
                else:
                    output_array[0][i] = self.NS_mass(MZAMS)
                    output_array[1][i] = 2
            elif (MZAMS > 18.5) and (MZAMS <= 21.7):
                if n > 50:
                    output_array[0][i] = self.BH_mass_low(MZAMS, 0.9)
                    output_array[1][i] = 3
                else:
                    output_array[0][i] = self.NS_mass(MZAMS)
                    output_array[1][i] = 2
            elif (MZAMS > 21.7) and (MZAMS <= 25.2):
                output_array[0][i] = self.BH_mass_low(MZAMS, 0.9)
                output_array[1][i] = 3
            elif (MZAMS > 25.2) and (MZAMS <= 27.5):
                if n > 65:
                    output_array[0][i] = self.BH_mass_low(MZAMS, 0.9)
                    output_array[1][i] = 3
                else: 
                    output_array[0][i] = self.NS_mass(MZAMS)
                    output_array[1][i] = 2
            elif (MZAMS > 27.5) and (MZAMS <= 60):
                if MZAMS > 42.22:
                    output_array[0][i] = self.BH_mass_high(MZAMS)
                    output_array[1][i] = 3
                else:
                    output_array[0][i] = self.BH_mass_low(MZAMS, 0.9)
                    output_array[1][i] = 3
            elif (MZAMS > 60) and (MZAMS <= 120):
                if n > 40:
                    output_array[0][i] = self.BH_mass_high(MZAMS)
                    output_array[1][i] = 3
                else:
                    output_array[0][i] = self.NS_mass(MZAMS)
                    output_array[1][i] = 2
            elif MZAMS < 0.5:
                output_array[0][i] = -99
                output_array[1][i] = -1
            else:
                output_array[0][i] = -99
                output_array[1][i] = -1
        return output_array
