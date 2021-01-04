import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
import os
import hoki
from hoki import load
import glob
import unittest

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
import os
import hoki
from hoki import load
import glob
possible_secondary_q = ['0.1', '0.2', '0.3',
                        '0.4', '0.5',
                        '0.6', '0.7', '0.8', '0.9']
# The list encompassing possible primary
# mass for primary and secondary. Trying to encompass
# all possible cases, even if there is a bit of redundancy and excess.
mass_list = [round(0.1, 1)] + [round(0.12 + x * 0.020, 2) for x in range(95)]
mass_list = mass_list + [round(2.05 + x * 0.05, 2) for x in range(20)]
mass_list = mass_list + [round(3.1 + 0.1 * x, 1) for x in range(70)]
mass_list = (mass_list + [11 + x for x in range(90)] +
             [120, 400, 500] + [125 + 25 * x for x in range(8)])


def convert(x):
    """
    Helper function to use for possible mass list. If the
    number can be represented by an integer,
    it turns into an integer.
    This is to make sure that I am really matching masses.
    """
    if (x % 1 == 0.0):
        return int(x)
    else:
        return x


mass_list = list(map(convert, mass_list))
period_list = ['0', '0.2', '0.4', '0.6', '0.8', '1', '1.2', '1.4',
               '1.6', '1.8', '2', '2.2', '2.4', '2.6', '2.8', '3',
               '3.2', '3.4', '3.6', '3.8', '4']
combos = []
for y in period_list:
    for z in possible_secondary_q:
        combos.append((z, y))


def assign_props(dictionary, input_str, y):
    dictionary[input_str] = y
    return input_str

# Did I end up scrambling the Naming Convention?
# Fitsmodels<sneplot-prefixed Filename><sin, bin, sec, or hmg>.fits
# May seem like an obvious one, but put it anyway:
# Helps ensure that we are not creating jiberrish file names.
class test_reformatter(unittest.TestCase):
    def help_test2(self, destination, metallicity):
        x = metallicity
        x1 = glob.glob("{}/{}/*sin.fits".format(destination, metallicity))
        x2 = glob.glob("{}/{}/*bin.fits".format(destination, metallicity))
        x3 = glob.glob("{}/{}/*hmg.fits".format(destination, metallicity))
        x4 = glob.glob("{}/{}/*sec.fits".format(destination, metallicity))
        total_set = set(x1+x2+x3+x4)
        g = glob.glob("{}/NEWSINMODS/{}/*".format(hoki.MODELS_PATH, x))
        h  = glob.glob("{}/NEWSINMODS/{}hmg/*".format(hoki.MODELS_PATH, x))
        j = glob.glob("{}/NEWBINMODS/NEWBINMODS/{}/*".format(hoki.MODELS_PATH, x))
        k = glob.glob("{}/NEWBINMODS/NEWSECMODS/{}_2/*".format(hoki.MODELS_PATH, x))
        for x in total_set:
            if (x[-8:]=="sin.fits"):
                if (not (os.path.isfile("{}/NEWSINMODS/{}/{}".format(hoki.MODELS_PATH, metallicity, x[len(destination + "/xxx/FitsModelss"):-8])))):
                    print("ERROR IN ", x)
                    return
            if (x[-8:]=="hmg.fits"):
                if not (os.path.isfile("{}/NEWSINMODS/{}hmg/{}".format(hoki.MODELS_PATH, metallicity, x[len(destination + "/xxx/FitsModelss"):-8]))):
                    print("Error in ", x)
                    return
            if (x[-8:]=="bin.fits"):
                if not (os.path.isfile("{}/NEWBINMODS/NEWBINMODS/{}/{}".format(hoki.MODELS_PATH, metallicity, x[len(destination + "/xxx/FitsModelss"):-8]))):
                    print("Error in ", x)
                    return
            if (x[-8:]=="sec.fits"):
                if not (os.path.isfile("{}/NEWBINMODS/NEWSECMODS/{}_2/{}".format(hoki.MODELS_PATH, metallicity, x[len(destination + "/xxx/FitsModelss"):-8]))):
                    print("Error in ", x)
                    return
        print(len(x1))
        print(len(g))
        print(len(x1)==len(g))
        print(len(x2))
        print(len(j))
        print(len(x2)==len(j))
        print(len(x3))
        print(len(h))
        print((len(x3)==len(h)))
        print(len(x4))
        print(len(k))
        print((len(x4)==len(k)))
        if (len(x1)==len(g)) and (len(x2)==len(j)) and (len(x3)==len(h)) and (len(x4)==len(k)):
            print("Passed! "+metallicity)
            return True
        else:
            print()
            print("Failed! "+metallicity)
            return False

    def extractor(self, metallicity, input_dir):
        """
            Portions of the extractor function meant to test a)
            whether all possible stars are covered
            by the looping hashmap-parameter matching scheme
            for binaries and single stars and b) whether all
            possible stars are getting covered by my assumptions
            about the length of 
        """
        # Set of files that exist and have been caught by the looping.
        caught_no = set()
        # Mapping of file name to a tuple of properties
        # (primary star mass,(secondary star mass, initial logP for bin.fits types)
        names_to_prop = {}
        # Find all filenames of reformatted BPASS models from reformatter
        for x in mass_list:
            # Find all NEWBINMODS systems of the specified metallicity
            for y in combos:
                st = "{}/{}/FitsModelssneplot-{}-{}-{}-{}bin.fits"
                st = st.format(input_dir,metallicity,
                               metallicity,
                               str(x), y[0],
                               y[1])
                names_to_prop[st] = (str(x), y[0], y[1])
                if (os.path.isfile(st)):
                    caught_no.add(st)
                # Find all single star systems of the specified metallicity
            st = "{}/{}/FitsModelssneplot-{}-{}sin.fits"
            st = st.format(input_dir, metallicity,
                               metallicity, str(x))
            names_to_prop[st] = (str(x), )
            if (os.path.isfile(st)):
                caught_no.add(st)
            # Find all Binary QHE systems of the specified metallicity
            st = "{}/{}/FitsModelssneplot-{}-{}hmg.fits"
            st = st.format(input_dir, metallicity,
                               metallicity, str(x))
            names_to_prop[st] = (str(x), )
            if (os.path.isfile(st)):
                caught_no.add(st)
            # Find all NEWSINMODS systems of the specified metallicity
            st = "{}/{}/FitsModelssneplot_2-{}-{}-*sec.fits"
            st = st.format(input_dir, metallicity,
                               metallicity,str(x))
            li = glob.glob(str(st))
            sec_files = set(li)
            caught_no = caught_no.union(sec_files)
            [assign_props(names_to_prop, name, (x,)) for name in sec_files]
        len_of_heading = len("{}/{}/FitsModelssneplot_2-{}-".format(input_dir, metallicity, metallicity))
        entries = glob.glob(str("{}/{}/*".format(input_dir, metallicity)))
        suffix_len = len("xxx.fits")  # 8
        # Note some of the branching statements have
        # if True statements.
        # This is redundant, but it would have been
        # a real pain to adjust the indentation
        # which in Python is really important
        # and quite annoying at times.
        # Count variable stands for number of sec.fits
        # files caught.
        count = 0
        for x in entries:
            rest = names_to_prop[x]
            if (x[-8:-5] == 'sec'):
                initlMass = float(rest[0])
                if (x[-10 - 1 * suffix_len] == "-"):
                    # Here the decimal is number is 9
                    # characters long and
                    # xxx.fits has length of 8
                    log_P_in_days = float(x[-9 - suffix_len:
                                                        -1 * suffix_len])
                    initlMass2 = float(x[(len_of_heading +
                                          len(str(rest[0]) +
                                          "-")): -10 - 1 * suffix_len])
                    st = "{}/{}/FitsModelssneplot_2-{}-{}-{}-{}sec.fits"
                    st = st.format(input_dir, metallicity,
                                   metallicity,rest[0], x[(len_of_heading +
                                                           len(str(rest[0]) +
                                                               "-")):
                                                          -10 - 1 * suffix_len],
                                   x[-9 - suffix_len: -1 * suffix_len])
                    if not ((log_P_in_days<0) and (os.path.isfile(st))):
                        return False
                    count+=1
                else:
                    # See if the number begins right
                    # after a dash in index -17
                    # 8 for the xxx.fits and 8 for
                    # the log_P and we want
                    # the preceding dash.
                    if (x[-9 - 1 * suffix_len] == "-"):
                        # Here the decimal is number is 8
                        # characters long and
                        # xxx.fits has length of 8.
                        log_P_in_days = float(x[-8 - suffix_len:
                                                -1 * suffix_len])
                        # Accounting for the dashes, find
                        # the sandwiched initial mass of the
                        # "remnant primary"
                        # sandwiched between initial secondary star
                        # mass and the initial logP
                        initlMass2 = float(x[(len_of_heading +
                                              len(str(rest[0]) +
                                                  "-")):
                                             -9 - 1 * suffix_len])
                        # Assume that the number would begin
                        # at index -16. (From my observation, we only
                        # have a few choices for number of digits of
                        # the log_period
                        st = "{}/{}/FitsModelssneplot_2-{}-{}-{}-{}sec.fits"
                        st = st.format(input_dir, metallicity,
                                       metallicity,rest[0], x[(len_of_heading +
                                                               len(str(rest[0]) +
                                                                   "-")):
                                                              -9 - 1 * suffix_len],
                                       x[-8 - suffix_len: -1 * suffix_len])
                        if not (os.path.isfile(st)):
                            return False
                        count+=1
                    elif (x[-8 - 1 * suffix_len] == "-"):
                        # Here the decimal is number is 7
                        # characters long and
                        # xxx.fits has length of 8
                        log_P_in_days = float(x[-7 - suffix_len:
                                                -1 * suffix_len])
                        # Accounting for the dashes on both
                        # ends of the "remnant primary" mass, find the
                        # sandwiched initial mass of the remnant primary
                        # sandwiched between initial secondary star mass
                        # and the initial logP
                        initlMass2 = float(x[(len_of_heading +
                                              len(str(rest[0]) +
                                                  "-")):
                                             -8 - 1 * suffix_len])
                        st = "{}/{}/FitsModelssneplot_2-{}-{}-{}-{}sec.fits"
                        st = st.format(input_dir, metallicity,
                                       metallicity,rest[0], x[(len_of_heading +
                                                               len(str(rest[0]) +
                                                                   "-")):
                                                              -8 - 1 * suffix_len],
                                       x[-7 - suffix_len: -1 * suffix_len])
                        if not (os.path.isfile(st) and (initlMass2>0)):
                            return False
                        count += 1
        return (count == len(glob.glob("{}/{}/FitsModelssneplot_2-{}-*sec.fits"
                                               .format(input_dir, metallicity,
                                                       metallicity))))
    def extractor_check_all(self, metallicity, input_dir):
        """
            Portions of the extractor function meant to test a)
            whether all possible stars are covered
            by the looping hashmap-parameter matching scheme
            for binaries and single stars and b) whether all
            possible stars are getting covered by my assumptions
            about the length of 
        """
        # Set of files that exist and have been caught by the looping.
        caught_no = set()
        # Mapping of file name to a tuple of properties
        # (primary star mass,(secondary star mass, initial logP for bin.fits types)
        names_to_prop = {}
        # Find all filenames of reformatted BPASS models from reformatter
        for x in mass_list:
            # Find all NEWBINMODS systems of the specified metallicity
            for y in combos:
                st = "{}/{}/FitsModelssneplot-{}-{}-{}-{}bin.fits"
                st = st.format(input_dir,metallicity,
                               metallicity,
                               str(x), y[0],
                               y[1])
                names_to_prop[st] = (str(x), y[0], y[1])
                if (os.path.isfile(st)):
                    caught_no.add(st)
                # Find all single star systems of the specified metallicity
            st = "{}/{}/FitsModelssneplot-{}-{}sin.fits"
            st = st.format(input_dir, metallicity,
                               metallicity, str(x))
            names_to_prop[st] = (str(x), )
            if (os.path.isfile(st)):
                caught_no.add(st)
            # Find all Binary QHE systems of the specified metallicity
            st = "{}/{}/FitsModelssneplot-{}-{}hmg.fits"
            st = st.format(input_dir, metallicity,
                               metallicity, str(x))
            names_to_prop[st] = (str(x), )
            if (os.path.isfile(st)):
                caught_no.add(st)
            # Find all NEWSINMODS systems of the specified metallicity
            st = "{}/{}/FitsModelssneplot_2-{}-{}-*sec.fits"
            st = st.format(input_dir, metallicity,
                               metallicity,str(x))
            li = glob.glob(str(st))
            sec_files = set(li)
            caught_no = caught_no.union(sec_files)
            [assign_props(names_to_prop, name, (x,)) for name in sec_files]
        return caught_no == set(glob.glob("{}/{}*.fits".format(input_dir, metallicity)))

    def test_reformatter_coverage(self, destination="/g/lu/scratch/ryotainagaki/BPASS_tester_newReformatTest/"):
        mets = ["zem5", "zem4", "zem3", "z001",
                "z002", "z003", "z004", "z006",
                "z008", "z010", "z014", "z020",
                "z030", "z040"]
        for met in mets:
            print(met)
            self.assertTrue(self.help_test2(destination, met))
    def test_extractor(self, source="/g/lu/scratch/ryotainagaki/BPASS_tester_newReformatTest/"):
        print("Testing if all of my assumptions regarding the reading of " +
               "NEWSECMODS type reformatted files are correct.")
        print("Also testing for coverage of ALL secmods files" +
              " in the sense that all of them would be considered" +
              " in isochrone file making.")
        mets = ["zem5", "zem4", "zem3", "z001",
                "z002", "z003", "z004", "z006",
                "z008", "z010", "z014", "z020",
                "z030", "z040"]
        for met in mets:
            self.assertTrue(self.extractor(met, source))
    def test_extractor2(self, source="/g/lu/scratch/ryotainagaki/BPASS_tester_newReformatTest/"):
        print("Also testing for whether all files in the input directory will be covered" +
              " by the infrastructure/assumptions of the extractor function")
        mets = ["zem5", "zem4", "zem3", "z001",
                "z002", "z003", "z004", "z006",
                "z008", "z010", "z014", "z020",
                "z030", "z040"]
        for met in mets:
            self.assertTrue(self.extractor_check_all(met, source))
if __name__ == '__main__':
    unittest.main()
