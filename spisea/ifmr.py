##########################################################
#
#
#IFMR_Raithel17 comes from Raithel et al. 2017 and has no metallicity dependence
#https://arxiv.org/pdf/1712.00021.pdf
#
#IFMR_Spera15 comes from Spera et al. 2015 appendix C and includes metallicity dependence
#https://arxiv.org/pdf/1505.05201.pdf
#
#########################################################

import numpy as np

class IFMR(object):
    def __init__(self):
        pass

    def get_Z(self, Fe_H):
        """
        This function converts metallicity given as [Fe/H] into Z values.
        It comes from Sharma et al. 2011 "Galaxia: A Code To Generate A Synthetic Survey of The Milky Way" equation (15) 
        https://ui.adsabs.harvard.edu/abs/2011ascl.soft01007S/abstract

        """
        return 10**(Fe_H - 1.72125)

    def Kalirai_mass(self, MZAMS):
        """                                                                                                      
        From Kalirai+07                                                                                          
        1.16 < MZAMS < 6.5 
        But we use this function for anything between 0.5 and 9 depending on the IFMR.
        FIXME: need to extend these ranges... explain extension somewhere? Paper maybe?                                                                   
        """

        result = 0.109*MZAMS + 0.394

        final = np.zeros(len(MZAMS))
        
        bad_idx = np.where(MZAMS < 0.5)
        final[bad_idx] = -99
        
        good_idx = np.where(MZAMS >= 0.5)
        final[good_idx] = result[good_idx]

        return final
        
        


class IFMR_Spera15(IFMR):
    """
    This IFMR comes from Spera et. al. 2015 Appendix C (Used for all MZAMS>= 7 M_sun)
    https://ui.adsabs.harvard.edu/abs/2015MNRAS.451.4086S/abstract

    The WD IFMR (used for MZAMS< 7_Msun) comes from 
    `Kalirai et al. (2008) <https://ui.adsabs.harvard.edu/abs/2008ApJ...676..594K/abstract>`_
    
    """
    
    #The get_Mco functions come from C11 of Spera
    def get_Mco_low_metal(self, Z, MZAMS):
        """
        C15 of Spera, valid for Z < 1.0e-3
    
        """
        
        B1 = 67.07
        K1 = 46.89
        K2 = 1.138e2
        d1 = 2.199e-2
        d2 = 2.602e-2


        g1 = 0.5/(1+10**((K1-MZAMS)*(d1))) #C12 of Spera
        g2 = 0.5/(1+10**((K2-MZAMS)*(d2))) #C12 of Spera

        return -2.0 + (B1 + 2)*(g1 + g2) #C11 of Spera


    def get_Mco_med_metal(self, Z, MZAMS):
        """
        C14 of Spera, valid for Z >= 1.0e-3 and Z <= 4.0e-3

        """
        
        B1 = 40.98 + 3.415e4*Z - 8.064e6*Z**2
        K1 = 35.17 + 1.548e4*Z - 3.759e6*Z**2
        K2 = 20.36 + 1.162e5*Z - 2.276e7*Z**2
        d1 = 2.500e-2 - 4.346*Z + 1.340e3*Z**2
        d2 = 1.750e-2 + 11.39*Z - 2.902e3*Z**2


        g1 = 0.5/(1+10**((K1-MZAMS)*(d1))) #C12 of Spera
        g2 = 0.5/(1+10**((K2-MZAMS)*(d2))) #C12 of Spera

        return -2.0 + (B1 + 2.0)*(g1 + g2) #C11 of Spera


    def get_Mco_high_metal(self, Z, MZAMS):
        """
        C13 of Spera, valid for Z > 4.0e-3

        """
        
        B1 = 59.63 - 2.969e3*Z + 4.988e4*Z**2
        K1 = 45.04 - 2.176e3*Z + 3.806e4*Z**2
        K2 = 1.389e2 - 4.664e3*Z + 5.106e4*Z**2
        d1 = 2.790e-2 - 1.780e-2*Z + 77.05*Z**2
        d2 = 6.730e-3 + 2.690*Z - 52.39*Z**2

        g1 = 0.5/(1+10**((K1-MZAMS)*(d1))) #C12 of Spera
        g2 = 0.5/(1+10**((K2-MZAMS)*(d2))) #C12 of Spera

        
        return -2.0 + (B1 + 2.0)*(g1 + g2) #C11 of Spera


    def get_Mco(self, Z, MZAMS):
        """
        This function uses Spera C11-C15 in order to reurn an array of core masses from an array of metallicities
        and ZAMS masses. It will be the same length as these two arrays with -99 entries where invalid (ie MZAMS< 7 M_sun)

        Parameters: 
        
        Z: an array with metallicities reported as Z where Z is metal_mass/total_mass
        MZAMS: an array of ZAMS masses in solar masses. The Spera functions are valid for MZAMS> 7 M_sun

        """
        #intialize an array of core masses with all entries equal to zero
        core_masses = np.zeros(len(MZAMS))

        #assign masses outside the range of validity for Spera a value of -99
        invalid_idx = np.where(MZAMS < 7.0)
        core_masses[invalid_idx] = -99

        #assign stars with Z < 1.0e-3 a core mass using the low metallicity core mass function get_Mco_low_metal
        low_metal_idx = np.where((Z < 1.0e-3) & (MZAMS >= 7.0))
        core_masses[low_metal_idx] = self.get_Mco_low_metal(Z[low_metal_idx], MZAMS[low_metal_idx])

        #assign stars with 1.0e-3 <= Z <= 4.0e-3 a core mass using the medium metallicity core mass function get_Mco_med_metal
        med_metal_idx = np.where((Z <= 4.0e-3) & (Z >= 1.0e-3) & (MZAMS >= 7.0))
        core_masses[med_metal_idx] = self.get_Mco_med_metal(Z[med_metal_idx], MZAMS[med_metal_idx])

        #assign stars with Z > 4.0e-3 a core mass using the high metallicity core mass function get_Mco_high_metal
        high_metal_idx = np.where((Z > 4.0e-3) & (MZAMS >= 7.0))
        core_masses[high_metal_idx] = self.get_Mco_high_metal(Z[high_metal_idx], MZAMS[high_metal_idx])

        return core_masses
        


    def M_rem_very_low_metal_low_mass(self, Z, Mco):
        """
        C1 of Spera, valid for Z <= 5.0e-4 and Mco <= 5.0
        
        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """
        
        p = -2.333 + 0.1559*Mco + 0.2700*Mco**2 #C2 of Spera

        #need to return p or 1.27, whichever is greater
        final = np.zeros(len(p))
        
        p_max_idx = np.where(p >= 1.27)
        final[p_max_idx] = p[p_max_idx]
        
        p_min_idx = np.where(p < 1.27)
        final[p_min_idx] = 1.27

        return final


    def M_rem_very_low_metal_med_mass(self, Z, Mco):
        """
        C1 of Spera, valid for Z <= 5.0e-4 and 5.0 < Mco < 10.0

        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """

        p = -2.333 + 0.1559*Mco + 0.2700*Mco**2 #C2 of Spera

        return p

    def M_rem_very_low_metal_high_mass(self, Z, Mco):
        """

        C1 of Spera, valid for Z <= 5.0e-4 and Mco >= 10.0

        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """

        p = -2.333 + 0.1559*Mco + 0.2700*Mco**2 #C2 of Spera

        m = -6.476e2*Z + 1.911 #C3 of Spera
        q = 2.300e3*Z + 11.67 #C3 of Spera
            
        f = m*Mco + q #C2 of Spera

        #need to return either p or f, whichever is less
        final = np.zeros(len(p))

        p_min_idx = np.where(p <= f)
        final[p_min_idx] = p[p_min_idx]

        f_min_idx = np.where(f < p)
        final[f_min_idx] = f[f_min_idx]
        
        return final

    def M_rem_low_metal_low_mass(self, Z, Mco):
        """
        C4 of Spera, valid for 5.0e-4 < Z < 1.0e-3, and Mco <= 5.0

        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """

        #values from C7 of Spera
        A1 = 1.105e5*Z - 1.258e2
        A2 = 91.56 - 1.957e4*Z - 1.558e7*Z**2
        L = 1.134e4*Z - 2.143
        n = 3.090e-2 - 22.30*Z + 7.363e4*Z**2

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        #need to return h or 1.27, whichever is greater
        final = np.zeros(len(h))
        
        h_max_idx = np.where(h >= 1.27)
        final[h_max_idx] = h[h_max_idx]
        
        h_min_idx = np.where(h < 1.27)
        final[h_min_idx] = 1.27

        return final
    

    def M_rem_low_metal_med_mass(self, Z, Mco):
        """
        C4 of Spera, valid for 5.0e-4 < Z < 1.0e-3, and 5.0 < Mco < 10.0

        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """

        #values from C7 of Spera
        A1 = 1.105e5*Z - 1.258e2
        A2 = 91.56 - 1.957e4*Z - 1.558e7*Z**2
        L = 1.134e4*Z - 2.143
        n = 3.090e-2 - 22.30*Z + 7.363e4*Z**2

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        return h

    def M_rem_low_metal_high_mass(self, Z, Mco):
        """
        C4 of Spera, valid for 5.0e-4 < Z < 1.0e-3, and Mco >= 10.0

        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """

        #values from C7 of Spera
        A1 = 1.105e5*Z - 1.258e2
        A2 = 91.56 - 1.957e4*Z - 1.558e7*Z**2
        L = 1.134e4*Z - 2.143
        n = 3.090e-2 - 22.30*Z + 7.363e4*Z**2

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        #values from C10 of Spera
        m = -6.476e2*Z + 1.911
        q = 2.300e3*Z + 11.67

        f = m*Mco + q #C5 of Spera

        #need to return either h or f, whichever is greater
        final = np.zeros(len(h))

        h_max_idx = np.where(h >= f)
        final[h_max_idx] = h[h_max_idx]

        f_max_idx = np.where(f > h)
        final[f_max_idx] = f[f_max_idx]
        
        return final
            

    def M_rem_med_metal_low_mass(self, Z, Mco):
        """
        C4 of Spera, valid for 1.0e-3 <= Z <= 4.0e-3, and Mco <= 5.0

        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609))

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        #need to return h or 1.27, whichever is greater
        final = np.zeros(len(h))
        
        h_max_idx = np.where(h >= 1.27)
        final[h_max_idx] = h[h_max_idx]
        
        h_min_idx = np.where(h < 1.27)
        final[h_min_idx] = 1.27

        return final

    def M_rem_med_metal_med_mass(self, Z, Mco):
        """
        C4 of Spera, valid for 1.0e-3 <= Z <= 4.0e-3, and 5.0 < Mco < 10.0

        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609))

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        return h

    def M_rem_med_metal_high_mass_1(self, Z, Mco):
        """
        C4 of Spera, valid for 1.0e-3 <= Z < 2.0e-3, and Mco >= 10.0

        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609))

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        #values from C9 of Spera
        m = -43.82*Z + 1.304
        q = -1.296e4*Z + 26.98

        f = m*Mco + q #C5 of Spera
            
        #need to return either h or f, whichever is greater
        final = np.zeros(len(h))

        h_max_idx = np.where(h >= f)
        final[h_max_idx] = h[h_max_idx]

        f_max_idx = np.where(f > h)
        final[f_max_idx] = f[f_max_idx]
        
        return final

    def M_rem_med_metal_high_mass_2(self, Z, Mco):
        """
        C4 of Spera, valid for 2.0e-3 <= Z <= 4.0e-3, and Mco >= 10.0

        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609))

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        #values from C8 of Spera
        m = 1.217
        q = 1.061

        f = m*Mco + q #C5 of Spera
            
        #need to return either h or f, whichever is greater
        final = np.zeros(len(h))

        h_max_idx = np.where(h >= f)
        final[h_max_idx] = h[h_max_idx]

        f_max_idx = np.where(f > h)
        final[f_max_idx] = f[f_max_idx]
        
        return final


    def M_rem_high_metal_low_mass(self, Z, Mco):
        """
        C4 of Spera, valid for Z > 4.0e-3, and Mco <= 5.0

        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609))

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        #need to return h or 1.27, whichever is greater
        final = np.zeros(len(h))
        
        h_max_idx = np.where(h >= 1.27)
        final[h_max_idx] = h[h_max_idx]
        
        h_min_idx = np.where(h < 1.27)
        final[h_min_idx] = 1.27

        return final

    def M_rem_high_metal_med_mass(self, Z, Mco):
        """
        C4 of Spera, valid for Z > 4.0e-3, and 5.0 < Mco < 10.0

        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609))

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        return h

    def M_rem_high_metal_high_mass(self, Z, Mco):
        """
        C4 of Spera, valid for Z > 4.0e-3, and Mco >= 10.0

        Parameters:

        Z: an array of metallicities reported as metal_mass/total_mass
        Mco: an arrray of core masses in M_sun

        """

        #values from C6 of Spera
        A1 = 1.340 - 29.46/(1 + (Z/(1.110e-3))**(2.361))
        A2 = 80.22 - 74.73*Z**(0.965)/(2.720e-3 + Z**(0.965))
        L = 5.683 + 3.533/(1 + (Z/(7.430e-3))**(1.993))
        n = 1.066 - 1.121/(1 + (Z/(2.558e-2))**(0.609))

        h = A1 + (A2 - A1)/(1 + 10**((L- Mco)*n)) #C5 of Spera

        #values from C8 of Spera
        m = 1.217
        q = 1.061

        f = m*Mco + q #C5 of Spera
            
        #need to return either h or f, whichever is greater
        final = np.zeros(len(h))

        h_max_idx = np.where(h >= f)
        final[h_max_idx] = h[h_max_idx]

        f_max_idx = np.where(f > h)
        final[f_max_idx] = f[f_max_idx]
        
        return final


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
        Range of validity: MZAMS > 0.5 
        
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

        #get core masses from MZAMS and metallicity
        core_mass= self.get_Mco(Z_array, mass_array)

        # where outside the validity of Spera on the low end use the Kaliari WD IFMR (ie where MZAMS < 7.0 M_sun)
        Kal_idx = np.where(core_mass < 0)
        rem_mass_array[Kal_idx] = self.Kalirai_mass(mass_array[Kal_idx])

        

        ##### very low metallicity Z < 5.0e-4
        
        #remnant masses of stars with Z < 5.0e-4 and Mco < 5.0
        very_low_metal_low_mass_idx = np.where((Z_array < 5.0e-4) & (core_mass < 5.0) & (core_mass >= 0))
        rem_mass_array[very_low_metal_low_mass_idx] = self.M_rem_very_low_metal_low_mass(Z_array[very_low_metal_low_mass_idx], core_mass[very_low_metal_low_mass_idx])

        #remnant masses of stars with Z < 5.0e-4 and 5.0 <= Mco <= 10.0
        very_low_metal_med_mass_idx = np.where((Z_array < 5.0e-4) & (core_mass >= 5.0) & (core_mass <= 10.0))
        rem_mass_array[very_low_metal_med_mass_idx] = self.M_rem_very_low_metal_med_mass(Z_array[very_low_metal_med_mass_idx], core_mass[very_low_metal_med_mass_idx])

        #remnant masses of stars with Z < 5.0e-4 and Mco > 10.0
        very_low_metal_high_mass_idx = np.where((Z_array < 5.0e-4) & (core_mass > 10.0))
        rem_mass_array[very_low_metal_high_mass_idx] = self.M_rem_very_low_metal_high_mass(Z_array[very_low_metal_high_mass_idx], core_mass[very_low_metal_high_mass_idx])

        #### low metallicity 5.0e-4 <= Z < 1.0e-3
        
        #remnant masses of stars with 5.0e-4 <= Z < 1.0e-3 and Mco < 5.0
        low_metal_low_mass_idx = np.where((Z_array >= 5.0e-4) & (Z_array < 1.0e-3) & (core_mass < 5.0) & (core_mass >= 0))
        rem_mass_array[low_metal_low_mass_idx] = self.M_rem_low_metal_low_mass(Z_array[low_metal_low_mass_idx], core_mass[low_metal_low_mass_idx])

        #remnant masses of stars with 5.0e-4 <= Z < 1.0e-3 and 5.0 <= Mco <= 10.0
        low_metal_med_mass_idx = np.where((Z_array >= 5.0e-4) & (Z_array < 1.0e-3) & (core_mass >= 5.0) & (core_mass <= 10.0))
        rem_mass_array[low_metal_med_mass_idx] = self.M_rem_low_metal_med_mass(Z_array[low_metal_med_mass_idx], core_mass[low_metal_med_mass_idx])

        #remnant masses of stars with 5.0e-4 <= Z < 1.0e-3 and Mco > 10.0
        low_metal_high_mass_idx = np.where((Z_array >= 5.0e-4) & (Z_array < 1.0e-3) & (core_mass > 10.0))
        rem_mass_array[low_metal_high_mass_idx] = self.M_rem_low_metal_high_mass(Z_array[low_metal_high_mass_idx], core_mass[low_metal_high_mass_idx])
                
        #### medium metallicity 1.0e-3 <= Z <= 4.0e-3
        
        #remnant masses of stars with 1.0e-3 <= Z <= 4.0e-3 and Mco < 5.0
        med_metal_low_mass_idx = np.where((Z_array >= 1.0e-3) & (Z_array <= 4.0e-3) & (core_mass < 5.0) & (core_mass >= 0))
        rem_mass_array[med_metal_low_mass_idx] = self.M_rem_med_metal_low_mass(Z_array[med_metal_low_mass_idx],core_mass[med_metal_low_mass_idx])

        #remnant masses of stars with 1.0e-3 <= Z <= 4.0e-3 and 5.0 <= Mco <= 10.0
        med_metal_med_mass_idx = np.where((Z_array >= 1.0e-3) & (Z_array <= 4.0e-3) & (core_mass >= 5.0) & (core_mass <= 10.0))
        rem_mass_array[med_metal_med_mass_idx] = self.M_rem_med_metal_med_mass(Z_array[med_metal_med_mass_idx], core_mass[med_metal_med_mass_idx])

        #remnant masses of stars with 1.0e-3 <= Z < 2.0e-3 and Mco > 10.0
        med_metal_high_mass_idx_1 = np.where((Z_array >= 1.0e-3) & (Z_array < 2.0e-3) & (core_mass > 10.0))
        rem_mass_array[med_metal_high_mass_idx_1] = self.M_rem_med_metal_high_mass_1(Z_array[med_metal_high_mass_idx_1], core_mass[med_metal_high_mass_idx_1])

        #remnant masses of stars with 2.0e-3 <= Z <= 4.0e-3 and Mco > 10.0
        med_metal_high_mass_idx_2 = np.where((Z_array >= 2.0e-3) & (Z_array <= 4.0e-3) & (core_mass > 10.0))
        rem_mass_array[med_metal_high_mass_idx_2] = self.M_rem_med_metal_high_mass_2(Z_array[med_metal_high_mass_idx_2], core_mass[med_metal_high_mass_idx_2])
        
        #### high metallicity Z > 4.0e-3
        
        #remnant masses of stars with Z > 4.0e-3 and Mco < 5.0
        high_metal_low_mass_idx = np.where((Z_array > 4.0e-3) & (core_mass < 5.0) & (core_mass >= 0))
        rem_mass_array[high_metal_low_mass_idx] = self.M_rem_high_metal_low_mass(Z_array[high_metal_low_mass_idx], core_mass[high_metal_low_mass_idx])

        #remnant masses of stars with Z > 4.0e-3 and 5.0 <= Mco <= 10.0
        high_metal_med_mass_idx = np.where((Z_array > 4.0e-3) & (core_mass >= 5.0) & (core_mass <= 10.0))
        rem_mass_array[high_metal_med_mass_idx] = self.M_rem_high_metal_med_mass(Z_array[high_metal_med_mass_idx], core_mass[high_metal_med_mass_idx])

        #remnant masses of stars with Z > 4.0e-3 and MZAMS > 10.0
        high_metal_high_mass_idx = np.where((Z_array > 4.0e-3) & (core_mass > 10.0))
        rem_mass_array[high_metal_high_mass_idx] = self.M_rem_high_metal_high_mass(Z_array[high_metal_high_mass_idx], core_mass[high_metal_high_mass_idx])
        
        #assign object types based on remnant mass
        bad_idx = np.where(rem_mass_array < 0) #outside the range of validity for the ifmr
        WD_idx = np.where((rem_mass_array <= 1.4) & (rem_mass_array >= 0 )) #based on the Chandresekhar limit
        NS_idx = np.where((rem_mass_array > 1.4) & (rem_mass_array <= 3.0)) #based on figures 15-17 of Spera 
        BH_idx = np.where(rem_mass_array > 3.0) #based on figures 15-17 of Spera

        output_array[0][bad_idx] = rem_mass_array[bad_idx]
        output_array[1][bad_idx] = -1
        
        output_array[0][WD_idx] = rem_mass_array[WD_idx]
        output_array[1][WD_idx] = codes['WD']

        output_array[0][NS_idx] = rem_mass_array[NS_idx]
        output_array[1][NS_idx] = codes['NS']

        output_array[0][BH_idx] = rem_mass_array[BH_idx]
        output_array[1][BH_idx] = codes['BH']

        return output_array
        


class IFMR_Raithel17(IFMR):
    """
    This IFMR comes from Raithel et al. 2017
    https://arxiv.org/pdf/1712.00021.pd
    The IFMR is a combination of the 
    WD IFMR from 
    `Kalirai et al. (2008) <https://ui.adsabs.harvard.edu/abs/2008ApJ...676..594K/abstract>`_
    and the NS/BH IFMR from
    `Raithel et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018ApJ...856...35R/abstract>`_.

    See Lam et al. (submitted) for more details.
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
        Drawing the mass from gaussian created using observational data
        FIXME: Write up in sam's paper? I (sam) got the distribution from Emily Ramey
        """
        return np.random.normal(loc=1.36, scale=0.09, size=len(MZAMS))

    def generate_death_mass(self, mass_array):
        """
        The top-level function that assigns the remnant type 
        and mass based on the stellar initial mass. 
        
        Parameters
        ----------
        mass_array: array of floats
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
        output_array[0][id_array1] = self.Kalirai_mass(mass_array[id_array1])
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
