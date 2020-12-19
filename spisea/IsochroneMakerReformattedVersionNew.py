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


# Source: https://stackoverflow.com/questions/44369504/
# how-to-convert-entire-dataframe-values-to-float-in-pandas
vals = hoki.dummy_dict.values()
vals = [hoki.dummy_dict[keyword] for keyword in ['timestep', 'age',
                                                 'log(R1)', 'log(T1)',
                                                 'log(L1)', 'M1', 'X',
                                                 'P_bin', 'log(a)', 'M2',
                                                 'log(R2)', 'log(T2)',
                                                 'log(L2)']]
cols_to_keep = ["col" + str(v + 1) for v in vals]
# According to the BPASS v2.2.1 manual, much of
# columns 50 and onward
# are basically spectra and atmosphere related models.
# They will be calculated later.
# I have only kept the column numbers corresponding to
# Source for inverse mapping: https://stackoverflow.com/questions/483666/
# reverse-invert-a-dictionary-mapping
# I create a mapping from column NUMBER to
# column name
invmap = {u: v for v, u in hoki.dummy_dict.items()}
# Accounting for zero indexing, write out a list
# of column names to assign to each BPASS dat file
# column name.
# Keep in mind the non-zero indexing of the data
# table itself.
lisnp = [invmap[int(x[3:]) - 1] for x in cols_to_keep]


def find_mergers(star_table, initial_mass):
    """
    given star_table, an Astropy Table of primary masses of a
    NEWBINMODS stellar model
    and age, and given the initial mass of the stellar model,
    sorts the star_table by age and
    returns the index at which the primary star's mass increases.
    If the mass of the primary star ends up being greater than the
    initial mass at the very first
    age represented by the star_table, then 0 is designated as
     the index where the merger occurs.
    If the merger does not occur then the length of the star_table
     is returned. (cannot be a valid index
    in a table). This is to ensure that all rows corresponding
     to positions after the returned index of the
    star_table are merger.
    """

    if (star_table['M1'][0] > initial_mass + 0.01):
        return 0
    # returns index where the merger is completed
    # assuming that f is a merger model.
    # We need to make sure that Mass of primary is increasing as age increases.
    # I make sure that age is actually increasing and
    # that the mass of the primary
    # is increasing by at least 0.01 solar masses.
    # Per Dr. Jan Edlridge, the increase in mass is supposed to be
    # an indicator of a merger happening
    val = np.where((np.diff(star_table['M1']) > 0.01) &
                   (np.diff(star_table['age']) >= 0))[0]
    if (len(val)):
        return val[-1] + 1
    return len(star_table)


def reformatter(destination, metallicity):
    """I will expedite the process of writing and saving files by
    making the BPASS stellar input files binary fits tables.
    This function reads every stellar model of the specified metallicities
    as an astropy table. Deletes photometry related columns of each stellar
    model and renames the columns to a human readable name as prescribed by
    the hoki.dummy_dict. Saves the reformatted BPASS stellar evolution model
    as a FITS formatted file for the sake of speed of processing later on.

    Parameters
    ----------
    destination: String
        describes the absolute path to the directory of BPASS models.
        Has a dash at the end
    metallicity: List of strings
        represents the metallicities of the BPASS input files that are
        to be processed.

    Remarks: Run this on Terminal not on Ipython. This function goes
    through tens of thousands of files full of data; expect it to take
    a long time to complete during your setup if you choose to use this
    via Jupyter Notebook. Using the reformatter on the the
    terminal will help save time.
        Cluster sub-class that produces a *resolved* stellar cluster.
    A table is output with the synthetic photometry and intrinsic
    properties of the individual stars (or stellar systems, if
    mutliplicity is used in the IMF object).

    If multiplicity is used, than a second table is produced that
    contains the properties of the companion stars independent of their
    primary stars.

    Parameters
    -----------
    iso: isochrone object
        SPISEA isochrone object

    imf: imf object
        SPISEA IMF object

    cluster_mass: float
        Total initial mass of the cluster, in M_sun

    ifmr: ifmr object or None
        If ifmr object is defined, will create compact remnants
        produced by the cluster at the given isochrone age. Otherwise,
        no compact remnants are produced.

    seed: int
        If set to non-None, all random sampling will be seeded with the
        specified seed, forcing identical output.
        Default None

    vebose: boolean
        True for verbose output.
    """
    # The list of possible BPASS metallicities.
    # HMG stands for Qasi-Chemically Homogeneous Evolution Models.
    # Real Deal metallicities=["zem5", "zem4", "z001",
    # "z002", "z003", "z004", "z006", "z008",
    # "z010", "z014","z020","z030", "z040", "zem5_hmg",
    # "zem4_hmg", "z001_hmg", "z002_hmg", "z003_hmg", "z004_hmg"]
    # I will use the Python HashSet since (assuming that strings
    # have a good HashCode, adding to the HashSet should be fast in
    # normal case.
    if (not os.path.isdir(destination)):
        os.makedirs(destination)
    setofAll = set()
    for x in metallicity:
        # Generate all possible BPASS models with the given metallicities.
        # Get a set of all single star models.
        a = set(glob.glob("{}/NEWSINMODS/{}/*".format(hoki.MODELS_PATH, x)))
        # Get a set of all Quasi-Chemicall Homogeneous Evolution (Binary)
        # Models.
        d = set(glob.glob("{}/NEWSINMODS/{}hmg/*".format(hoki.MODELS_PATH, x)))
        # Get a set of all regular binary and merged binary models.
        b = set(glob.glob(("{}/NEWBINMODS/" +
                           "NEWBINMODS/{}/*").format(hoki.MODELS_PATH, x)))
        # Get a set of all models where the primary is a compact remnant.
        c = set(glob.glob(("{}/NEWBINMODS/" +
                          "NEWSECMODS/{}_2/*").format(hoki.MODELS_PATH, x)))
        print(("{}/NEWBINMODS/"+"NEWSECMODS/{}_2/*").format(hoki.MODELS_PATH, x))
        print(len(c))
        setofAll = setofAll.union(a, b, c, d)
        if not os.path.isdir("{}/{}/".format(destination, x)):
            os.makedirs("{}/{}/".format(destination, x))
        # Creates container directory for stellar models grouped by metallicity
        # and whether they are single or binary related
        # the latter included mergers, QHE Models
        # Single star models that are not the hmg models
        # go into the <Metallicity>/sin
        # directory.
    new_sin_to_met = len(hoki.MODELS_PATH) + len("/NEWSINMODS/")
    new_sin_and_met = len(hoki.MODELS_PATH) + len("/NEWSINMODS/zxxx/")
    new_hmg_and_met = len(hoki.MODELS_PATH) + len("/NEWSINMODS/zxxxhmg/")
    new_bin_to_met = len(hoki.MODELS_PATH) + len("/NEBINMODS/NEWBINMODS/")
    new_sec_to_met = len(hoki.MODELS_PATH) + len("/NEBINMODS/NEWSECMODS/")
    new_bin_and_met = len(hoki.MODELS_PATH) + \
    len("/NEBINMODS/NEWBINMODS/zxxx/")
    new_sec_and_met = len(hoki.MODELS_PATH) + \
    len("/NEBINMODS/NEWBINMODS/zxxx_2/")
    for x in setofAll:
        astroTable = Table.read(x, format='ascii')
        # Drop all of the columns that are not defined by the
        # HOKI package provided mapping of column name to
        # column number.
        # We will be recalculating the photometry related stuff
        # Using zero indexing for list of columns,
        # All columns (except column number 12, 13, 36)
        # up to 48 inclusive will be kept
        # when we create the isochrone.
        astroTable.keep_columns(cols_to_keep)
        astroTable.rename_columns(astroTable.colnames, lisnp)
        astroTable.sort("age")
        # Will remove photometry related columns
        # Merger Related indicates whether the star is
        # a model with a merger in it.
        # Checking if the model is a single star model.
        sin_segment = x[new_sin_to_met:new_sin_to_met + 4]
        # SEE IF THE FILE IS IN THE NEWSINMODS directory
        if (x[len(hoki.MODELS_PATH) + 1:len(hoki.MODELS_PATH) +
                len("NEWSINMODS") + 1] == "NEWSINMODS" and x[new_sin_to_met +
                                                         4:new_sin_to_met +
                                                         7] != "hmg"):
            print(destination +
                             "/" + sin_segment + "/Fits" +
                             "Models{}sin.fits".format(
                                    x[new_sin_and_met:]))
            astroTable.write(destination +
                             "/" + sin_segment + "/Fits" +
                             "Models{}sin.fits".format(
                                    x[new_sin_and_met:]),
                             format='fits', overwrite=True)
        else:
            # Checking if the model is a HMG Binary model
            # (Still for some reason in SINMODS)
            if x[len(hoki.MODELS_PATH) + 1:len(hoki.MODELS_PATH) + 1 +
                 len("NEWSINMODS")] == "NEWSINMODS":
                print(destination + "/" +
                                 sin_segment + "/Fits" +
                                 "Models{}hmg.fits".format(x[new_hmg_and_met:]))
                astroTable.write(destination + "/" +
                                 sin_segment + "/Fits" +
                                 "Models{}hmg.fits".format(x[new_hmg_and_met:]),
                                 format='fits', overwrite=True)
                continue
                # Models that count as secondary star with compact remnants
                # go into the <metallicity>/bin subdirectory of
                # the destination.
            if x[len(hoki.MODELS_PATH) + 1:len(hoki.MODELS_PATH) + 1 +
                 len("NEWBINMODS/XXXXXXXXXX")] == "NEWBINMODS/NEWSECMODS":
                print(destination +
                      "/" + x[new_sec_to_met:
                              (new_sec_and_met) - 2] +
                      "/FitsModels{}sec.fits".format(x[new_sec_and_met + 1:]))
                astroTable.write(destination +
                                 "/" + x[new_sec_to_met:
                                         (new_sec_and_met) - 2] +
                                 "/FitsModels{}sec.fits".format(x[new_sec_and_met + 1:]),
                                 format='fits', overwrite=True)
            else:
                # Models that count as secondary star with compact remnants
                # go into the <metallicity>/ subdirectory of the destinatuion.
                print((destination + '/' +
                       x[new_bin_to_met: new_bin_and_met] +
                       "/FitsModels{}bin.fits".
                       format(x[new_bin_and_met + 1:])))
                astroTable.write(destination + '/' +
                                 x[new_bin_to_met: new_bin_and_met] +
                                 "/FitsModels{}bin.fits".format(x[new_bin_and_met + 1:]),
                                 format='fits', overwrite=True)


def extractor(age, metallicity, input_dir, bpass_evo_dir,
              margin):
    """
    The following extracts the preprocessed BPASS
    isochrone file
    from the pool of listed files. From the <metallicity> subdirectory
    of the input directory, the function inspects all of the files
    and sees if there are rows in the files that correspond to
    the specified age (log age to be precise).
    Checks for merger related (models_type=0) models
    and mark them using the functions merger_related and
    find_mergers respectively.



    parameters
    ----------
    age: float or double
          The age in years of the isochrone desired
    metallicity: string
          The BPASS metallicity of the isochrone desired
    input_dir: string
          The absolute path to the directory with the reformatted BPASS stellar
          models. The path should be to a folder whose name was recently
          passed into reformatter.
    bpass_evo_dir: string
          The absolute path to the directory where the outputs of this program,
          preprocessed isochrone files for BPASS, will
          be stored in subdirectories
          named after their metallicities.
    margin: float or double
             Contains the margin of error of base-10 logAge that will be
             tolerated for rows of stellar model that are included in
             the isochrone.

    Output:
    Saves into the specified destination
    directory a reformatted, FITS format version
    of a preprocessed BPASS isochrone of the specified age
    and metallicity whose individual
    rows correspond to stars of age closest to
    the specified age and within the specified
    margin of error of log(age of the star in years)

    Parameters added to the table:
    ‘mass’ - The column corresponding to the initial mass
             of the primary/single star.
    ‘mass2’ -  The column corresponding to the initial mass of
               the secondary. Set to nan if the system is not a binary
    ‘log_P’ - The column corresponding to the log(Period in days) of
              the binary. Set to nan if the system is not a binary
    ‘Single’ - The column that tells whether the stellar system is single
    ‘Merged?’ - Whether the binary star has been merged.

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
                       metallicity,
                       str(x))
        li = glob.glob(str(st))
        sec_files = set(li)
        caught_no = caught_no.union(sec_files)
        [assign_props(names_to_prop, name, (x,)) for name in sec_files]
    len_of_heading = len("{}/{}/FitsModelssneplot_2-{}-".format(input_dir, metallicity, metallicity))
    bigOne = None
    initial = True
    initlMass = np.nan
    initlMass2 = np.nan
    indicesOfInterest = None
    entries = glob.glob(str("{}/{}/*".format(input_dir, metallicity)))
    suffix_len = len("xxx.fits")  # 8
    for x in entries:
        # If I want to see the name of the file, uncomment below.
        # print(x)
        # x is the name of the reformatted stellar evolution file.
        # Rest carries the initial mass of the system and for NEWBINMODS
        # systems carries the secondary and log_P in days
        rest = names_to_prop[x]
        # The following conditional would have been that x is a path to a file,
        # but this seems to be always true.
        # as I am getting the entries from globs and
        # prior checks as to whether they are existing files.
        # We only want to consider the
        org = Table.read(x, format='fits')
        indicesOfInterest = np.where(np.abs(np.log10(org['age']) -
                                            age) <= margin)[0]
        f = org[indicesOfInterest]  # f stands for frame in DataFrame
        indicesOfInterest = np.array(indicesOfInterest)
        # If no stars have a log10(age) within margin of
        # the given lage in log-10 years, we
        # must skip over to the next model.
        if (len(f) != 0):
            # Find the star with age closest to the input log(Age of the star)
            filterDown = np.where(np.abs(f['age'] - 10 ** age) ==
                                  np.min(np.abs(f['age'] - 10 ** age)))[0]
            f = f[filterDown]
            indicesOfInterest = np.array([indicesOfInterest[filterDown][0]])
            if len(f) != 0:
                f = f[[0]]
                f['source'] = 0
                indexlen = 1
                if (x[-8:-5] == 'bin'):
                    f['single'] = np.repeat(False, indexlen)
                    initlMass = float(rest[0])
                    f['mass'] = np.repeat(initlMass, indexlen)
                    f['mass2'] = np.repeat(initlMass * float(rest[1]),
                                           indexlen)
                    f['initl_logP'] = np.repeat(float(rest[2]), indexlen)
                    # Now, for binaries, I check whether a model is
                    # a merger model
                    # is the same for all rows of the model
                    merge_pt = find_mergers(org[['age', 'M1']], initlMass)
                    # Find whether we should treat our model as
                    # one big star or still two stars
                    # will still keep initial parameters
                    # for the sake of working
                    # with the Duchene-Krauss distributions
                    f['mergered?'] = indicesOfInterest >= merge_pt
                    f['source'] = 1
                elif (x[-8:-5] == 'sec'):
                    # Recall that the ending of the file
                    # is going to be XXX.fits
                    # The XXX can be sec, hmg, bin, sin.
                    # I will explot the pattern that the log_P is
                    # going to be either 9 characters,
                    # 8 characters, or 7 characters
                    # long.
                    f['single'] = np.repeat(False, indexlen)
                    initlMass = float(rest[0])
                    # See if the number would begin
                    # right after a dash in index -18.
                    # 8 characters given for the xxx.fits
                    # 9 spaces for the decimal number.
                    # We want the dash before that as our cue.
                    if (x[-10 - 1 * suffix_len] == "-"):
                        # Here the decimal is number is 9
                        # characters long and
                        # xxx.fits has length of 8
                        log_P_in_days = float(x[-9 - suffix_len:
                                                -1 * suffix_len])
                        # Accounting for the dashes, find the
                        # sandwiched initial
                        # mass of the "remnant primary"
                        # sandwiched between initial secondary star mass
                        # and the initial logP
                        initlMass2 = float(x[(len_of_heading +
                                              len(str(rest[0]) + "-")):
                                             -10 - 1 * suffix_len])
                        f['source'] = 2
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
                            f['source'] = 3
                        else:
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
                            f['source'] = 4
                    f['mass'] = np.repeat(initlMass2, indexlen)
                    f['mass2'] = np.repeat(initlMass, indexlen)
                    temp1 = f['M2']
                    temp2 = f['log(T2)']
                    temp3 = f['log(R2)']
                    temp4 = f['log(L2)']
                    f['M2'] = f['M1']
                    f['log(R2)'] = f['log(R1)']
                    f['log(T2)'] = f['log(T1)']
                    f['log(L2)'] = f['log(L1)']
                    f['M1'] = temp1
                    f['log(R1)'] = temp3
                    f['log(T1)'] = temp2
                    f['log(L1)'] = temp4
                    f['initl_logP'] = np.repeat(log_P_in_days, indexlen)
                    f['mergered?'] = np.repeat(False, indexlen)
                    f['single'] = np.repeat(False, indexlen)
                elif (x[-1 * suffix_len: -5] == 'hmg'):
                    f['single'] = np.repeat(True, indexlen)
                    initlMass = float(rest[0])
                    # To be consistent with obtaining initial mass
                    # (or whichever value is closest) for the HMG Model).
                    # Per Eldridge the QHE models are single star models.
                    initlMass2 = np.nan
                    P_in_days = np.nan
                    f['mass'] = np.repeat(initlMass, indexlen)
                    f['mass2'] = np.repeat(initlMass2, indexlen)
                    f['initl_logP'] = np.repeat(P_in_days, indexlen)
                    f['mergered?'] = np.repeat(False, indexlen)
                    f['source'] = 5
                else:
                    f['single'] = np.repeat(True, indexlen)
                    # Single stars do not have companions.
                    initlMass = rest[0]
                    initlMass2 = np.nan
                    P_in_days = np.nan
                    f['mass'] = np.repeat(initlMass, indexlen)
                    f['mass2'] = np.repeat(initlMass2, indexlen)
                    f['initl_logP'] = np.repeat(np.nan, indexlen)
                    f['mergered?'] = np.repeat(False, indexlen)
                    f['source'] = 6
                if initial:
                    initial = False
                    bigOne = f.to_pandas()
                else:
                    bigOne = pd.concat([f.to_pandas(), bigOne])
    if not isinstance(bigOne, type(None)) and not (bigOne['age'].empty):
        bigOne = bigOne.apply(pd.to_numeric, errors='coerce')
        reduced = Table.from_pandas(bigOne)
        if not os.path.isdir('{}/iso/{}/'.format(bpass_evo_dir, metallicity)):
            os.makedirs('{}/iso/{}/'.format(bpass_evo_dir, metallicity))
        reduced.write("{}iso/{}/iso{}.fits".format(bpass_evo_dir,
                                                   metallicity, str(age)),
                      format='fits',
                      overwrite=True)
