class IFMR_Spera15(IFMR):

    #The get_Mco functions come from C11 of Spera
    def get_Mco_low_metal(self, Z, MZAMS): #C15 of Spera, valid for Z<1.0e-3

        B1 = 67.07
        K1 = 46.89
        K2 = 1.138e2
        d1 = 2.199e-2
        d2 = 2.602e-2


        g1 = 0.5/(1+10**((K1-MZAMS)*(d1))) #C12 of Spera
        g2 = 0.5/(1+10**((K2-MZAMS)*(d2))) #C12 of Spera

        return -2.0 + (B1 + 2)*g1 + g2 #C11 of Spera


    def get_Mco_med_metal(self, Z, MZAMS): #C14 of Spera, valid for Z >= 1.0e-3 and Z<=4.0e-3

        B1 = 40.98 + 3.415e4*Z - 8.064e6*Z**2
        K1 = 35.17 + 1.548e4*Z - 3.759e6*Z**2
        K2 = 20.36 + 1.162e5*Z - 2.276e7*Z**2
        d1 = 2.500e-2 - 4.346*Z + 1.340e3*Z**2
        d2 = 1.750e-2 + 11.39*Z - 2.902e3*Z**2


        g1 = 0.5/(1+10**((K1-MZAMS)*(d1))) #C12 of Spera
        g2 = 0.5/(1+10**((K2-MZAMS)*(d2))) #C12 of Spera

        return -2.0 + (B1 + 2.0)*g1 + g2 #C11 of Spera


    def get_Mco_high_metal(self, Z, MZAMS): #C13 of Spera, valid for Z > 4.0e-3

        B1 = 59.63 - 2.969e3*Z + 4.988e4*Z**2
        K1 = 45.04 - 2.176e3*Z + 3.806e4*Z**2
        K2 = 1.389e2 - 4.664e3*Z + 5.106e4*Z**2
        d1 = 2.790e-2 - 1.780e-2*Z + 77.05*Z**2
        d2 = 6.730e-3 + 2.690*Z - 52.39*Z**2

        g1 = 0.5/(1+10**((K1-MZAMS)*(d1))) #C12 of Spera
        g2 = 0.5/(1+10**((K2-MZAMS)*(d2))) #C12 of Spera

        
        return -2.0 + (B1 + 2)*g1 + g2 #C11 of Spera


    def M_rem_very_low_metal_low_mass(self, Z, MZAMS): #C1 of Spera, valid for Z <= 5.0e-4 and MZAMS <= 5.0

        Mco = self.get_Mco_low_metal(Z, MZAMS)
        p = -2.333 + 0.1559*Mco + 0.2700*Mco*Mco #C2 of Spera

        return max(p, 1.27)

    def M_rem_very_low_metal_med_mass(self, Z, MZAMS): #C1 of Spera, valid for Z <= 5.0e-4 and 5.0 < MZAMS < 10.0

        Mco = self.get_Mco_low_metal(Z, MZAMS)
        p = -2.333 + 0.1559*Mco + 0.2700*Mco*Mco #C2 of Spera

        return p

    def M_rem_very_low_metal_high_mass(self, Z, MZAMS): #C1 of Spera, valid for Z <= 5.0e-4 and MZAMS >= 10.0

        Mco = self.get_Mco_low_metal(Z, MZAMS)

        p = -2.333 + 0.1559*Mco + 0.2700*Mco*Mco #C2 of Spera

        m = -6.476e2*Z + 1.911 #C3 of Spera
        q = 2.300e3*Z + 11.67 #C3 of Spera
            
        f = m*Mco + q #C2 of Spera
        return min(p, f)

    def M_rem_low_metal_low_mass(self, Z, MZAMS): #C4 of Spera, valid for 5.0e-4 < Z < 1.0e-3, and MZAMS <= 5.0

        Mco = self.get_Mco_low_metal(Z, MZAMS)

        #values from C7 of Spera
        A1 = 1.105*e5*Z - 1.258e2
        A2 = 91.56 - 1.957e4*Z - 1.558e7*Z**2
        L = 1.134e4*Z - 2.143
        n = 3.090e-2 - 22.30*Z + 7.363e4*Z**2

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        return max(h, 1.27)

    def M_rem_low_metal_med_mass(self, Z, MZAMS): #C4 of Spera, valid for 5.0e-4 < Z < 1.0e-3, and 5.0 < MZAMS < 10.0

        Mco = self.get_Mco_low_metal(Z, MZAMS)

        #values from C7 of Spera
        A1 = 1.105*e5*Z - 1.258e2
        A2 = 91.56 - 1.957e4*Z - 1.558e7*Z**2
        L = 1.134e4*Z - 2.143
        n = 3.090e-2 - 22.30*Z + 7.363e4*Z**2

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        return h

    def M_rem_low_metal_high_mass(self, Z, MZAMS): #C4 of Spera, valid for 5.0e-4 < Z < 1.0e-3, and MZAMS >= 10.0

        Mco = self.get_Mco_low_metal(Z, MZAMS)

        #values from C7 of Spera
        A1 = 1.105*e5*Z - 1.258e2
        A2 = 91.56 - 1.957e4*Z - 1.558e7*Z**2
        L = 1.134e4*Z - 2.143
        n = 3.090e-2 - 22.30*Z + 7.363e4*Z**2

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        #values from C10 of Spera
        m = -6.476e2*Z + 1.911
        q = 2.300e3*Z + 11.67

        f = m*Mco + q #C5 of Spera
            
        return max(h, f)

    def M_rem_med_metal_low_mass(self, Z, MZAMS): #C4 of Spera, valid for 1.0e-3 <= Z <= 4.0e-3, and MZAMS <= 5.0

        Mco = self.get_Mco_med_metal(Z, MZAMS)

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609))

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        return max(h, 1.27)

    def M_rem_med_metal_med_mass(self, Z, MZAMS): #C4 of Spera, valid for 1.0e-3 <= Z <= 4.0e-3, and 5.0 < MZAMS < 10.0

        Mco = self.get_Mco_med_metal(Z, MZAMS)

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609))

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        return h

    def M_rem_med_metal_high_mass_1(self, Z, MZAMS): #C4 of Spera, valid for 1.0e-4 <= Z < 2.0e-3, and MZAMS >= 10.0

        Mco = self.get_Mco_med_metal(Z, MZAMS)

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609)

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        #values from C9 of Spera
        m = -43.82*Z + 1.304
        q = -1.296e4*Z + 26.98

        f = m*Mco + q #C5 of Spera
            
        return max(h, f)

    def M_rem_med_metal_high_mass_2(self, Z, MZAMS): #C4 of Spera, valid for 2.0e-3 <= Z <= 4.0e-3, and MZAMS >= 10.0

        Mco = self.get_Mco_med_metal(Z, MZAMS)

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609)

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        #values from C8 of Spera
        m = 1.217
        q = 1.061

        f = m*Mco + q #C5 of Spera
            
        return max(h, f)


    def M_rem_high_metal_low_mass(self, Z, MZAMS): #C4 of Spera, valid for Z > 4.0e-3, and MZAMS <= 5.0

        Mco = self.get_Mco_high_metal(Z, MZAMS)

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609))

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        return max(h, 1.27)

    def M_rem_high_metal_med_mass(self, Z, MZAMS): #C4 of Spera, valid for Z > 4.0e-3, and 5.0 < MZAMS < 10.0

        Mco = self.get_Mco_high_metal(Z, MZAMS)

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609))

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        return h

    def M_rem_med_metal_high_mass(self, Z, MZAMS): #C4 of Spera, valid for Z > 4.0e-3, and MZAMS >= 10.0

        Mco = self.get_Mco_med_metal(Z, MZAMS)

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609)

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        #values from C8 of Spera
        m = 1.217
        q = 1.061

        f = m*Mco + q #C5 of Spera
            
        return max(h, f)


        def generate_death_mass(self, mass_array, metallicity_array):
        """
        The top-level function that assigns the remnant type 
        and mass based on the stellar initial mass. 
        
        Parameters
        ----------
        mass_array: array of floats
            Array of initial stellar masses. Units are
            M_sun.

        metallicity_array: array of floats
            Array of metallicities in terms of [Fe/H]

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

        Range of validity: 0.5 < MZAMS < 120 #FIXME need to find out what the limits of Spera's IFMR are and assign values appropriately

        Returns
        -------
        output_arr: 2-element array
            output_array[0] contains the remnant mass, and 
            output_array[1] contains the typecode
        """

        #output_array[0] holds the remnant mass
        #output_array[1] holds the remnant type
        output_array = np.zeros((2, len(mass_array)))

        codes = {'WD': 101, 'NS': 102, 'BH': 103}

        #create array to store the remnant masses generated by Spera equations
        rem_mass_array = np.zeros(len(mass_array))

        #convert metallicity_array into a Z_array-changing from [Fe/H] to Z
        Z_array = np.zeros((len(metallicity_array)))
        metal_idx = np.where(metallicity_array != np.nan)
        Z_array[metal_idx] = self.get_Z(metallicity_array[metal_idx])

        ##### very low metallicity Z < 5.0e-4
        
        #remnant masses of stars with Z < 5.0e-4 and MZAMS < 5.0
        very_low_metal_low_mass_idx = np.where(Z_array < 5.0e-4 and mass_array < 5.0)
        rem_mass_array[very_low_metal_low_mass_idx] = self.M_rem_very_low_metal_low_mass(Z_array[very_low_metal_low_mass_idx], mass_array[very_low_metal_low_mass_idx])

        #remnant masses of stars with Z < 5.0e-4 and 5.0 <= MZAMS <= 10.0
        very_low_metal_med_mass_idx = np.where(Z_array < 5.0e-4 and mass_array >= 5.0 and mass_array <= 10.0)
        rem_mass_array[very_low_metal_med_mass_idx] = self.M_rem_very_low_metal_med_mass(Z_array[very_low_metal_med_mass_idx], mass_array[very_low_metal_med_mass_idx])

        #remnant masses of stars with Z < 5.0e-4 and MZAMS > 10.0
        very_low_metal_high_mass_idx = np.where(Z_array < 5.0e-4 and mass_array > 10.0)
        rem_mass_array[very_low_metal_high_mass_idx] = self.M_rem_very_low_metal_high_mass(Z_array[very_low_metal_high_mass_idx], mass_array[very_low_metal_high_mass_idx])

        ####low metallicity 5.0e-4 <= Z < 1.0e-3
        
        #remnant masses of stars with 5.0e-4 <= Z < 1.0e-3 and MZAMS < 5.0
        low_metal_low_mass_idx = np.where(Z_array >= 5.0e-4 and Z_array < 1.0e-3 and mass_array < 5.0)
        rem_mass_array[low_metal_low_mass_idx] = self.M_rem_low_metal_low_mass(Z_array[low_metal_low_mass_idx], mass_array[low_metal_low_mass_idx])

        #remnant masses of stars with 5.0e-4 <= Z < 1.0e-3 and 5.0 <= MZAMS <= 10.0
        low_metal_med_mass_idx = np.where(Z_array >= 5.0e-4 and Z_array < 1.0e-3 and mass_array >= 5.0 and mass_array <= 10.0)
        rem_mass_array[low_metal_med_mass_idx] = self.M_rem_low_metal_med_mass(Z_array[low_metal_med_mass_idx], mass_array[low_metal_med_mass_idx])

        #remnant masses of stars with 5.0e-4 <= Z < 1.0e-3 and MZAMS > 10.0
        low_metal_high_mass_idx = np.where(Z_array >= 5.0e-4 and Z_array < 1.0e-3 and mass_array > 10.0)
        rem_mass_array[low_metal_high_mass_idx] = self.M_rem_low_metal_high_mass(Z_array[low_metal_high_mass_idx], mass_array[low_metal_high_mass_idx])

        
        #assign object types based on remnant mass #FIXME-not sure I understand what Spera is saying. I feel like there shoud be some probablistic way of doing this instead of having hard cutoffs
        WD_idx = np.where(rem_mass_array <= 1.4 ) #based on the Chandresekhar limit
        NS_idx = np.where(rem_mass_array > 1.4 and rem_mass_array <= 3.) #based on figures 15-17 of Spera (dubious because neutron stars can be less than 1.4 Msun (1.25 Msun)
        BH_idx = np.where(rem_mass_array > 3.) #based on figures 15-17 of Spera (also dubious)

        output_array[0][WD_idx] = rem_mass_array[WD_IDX]
        output_array[1][WD_idx] = codes['WD']

        output_array[0][NS_idx] = rem_mass_array[NS_IDX]
        output_array[1][NS_idx] = codes['NS']

        output_array[0][BH_idx] = rem_mass_array[BH_IDX]
        output_array[1][BH_idx] = codes['BH']

        return output_array
