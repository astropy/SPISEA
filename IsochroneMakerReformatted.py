import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
import os
import hoki
from hoki import load
import glob
import re

"""Source: https://stackoverflow.com/questions/44369504/
how-to-convert-entire-dataframe-values-to-float-in-pandas"""

colsToElim = [x for x in range(96) if x not in [hoki.dummy_dict[y]
              for y in hoki.dummy_dict.keys()]]
def merger_related(f):
    # Code that eventually returns whether
    # The given model (in form of dataframe)
    # has the conditions necessary
    # to be a merger model.
    return False
    
def find_mergers(f):
    # returns when the merger occurs (first time)
    # assuming that f is a merger model.
    return False


def reformatter(bevo_dir):
    """I will expedite the process of writing and saving files by
    making the BPASS stellar input files binary fits tables This method
    creates a reformatted FITS version of 
    all of the available BPASS data so that data can be loaded quickly
    and processed in the extractor method."""

    # The list of possible BPASS metallicities.
    # HMG stands for Qasi-Chemically Homogeneous Evolution Models.

    metallicities=["zem5", "zem4", "z001", "z002", "z003", "z004", "z006", "z008",
                   "z010", "z014","z020","z030", "z040", "zem5_hmg",
                   "zem4_hmg", "z001_hmg", "z002_hmg", "z003_hmg", "z004_hmg"]
    # I will use the Python HashSet since (assuming that strings
    # have a good HashCode, adding to the HashSet should be fast

    setofAll = set()
    for x in metallicity:
        a = set(glob.glob(hoki.MODELS_PATH + "/NEWSINMODS/" + x + "/*"))
        b = set(glob.glob(hoki.MODELS_PATH + "/NEWBINMODS/NEWBINMODS/" + x + "/*"))
        c = set(glob.glob(hoki.MODELS_PATH + "/NEWBINMODS/NEWSECMODS/" + x + "_2/*"))
        setofAll = setofAll.union(a, b, c)
    
    for x in setofAll:

        f = pd.read_csv(hoki.MODELS_PATH + x, sep = "\s+",
                                header = None)
        f.apply(pd.to_numeric, errors = 'coerce')
        f.transpose()
        
        # Drop all of the columns that are not defined by the
        # HOKI package provided mapping of column name to
        # column number.
        f.drop(columns = [f.columns[x] for x in colsToElim],
               inplace = True)
        # Source: https://stackoverflow.com/questions/483666/
        # reverse-invert-a-dictionary-mapping
        # I create a mapping from column NUMBER to 
        # column name
        invmap = {u: v for v, u in hoki.dummy_dict.items()}
        f.rename(columns = invmap, inplace = True)
        f["model_type"] = [typ] * len(f.index)
        astroTable = Table.from_pandas(f)
        # We will be recalculating the photometry related stuff
        # when we create the isochrone.
        # Will remove photometry related columns
        astroTable.remove_columns(['V-I', 'U',
                                   'B', 'V', 'R', 'I', 'J', 'H', 'K','u',
                                   'g', 'r', 'i', 'z', 'f300w', 'f336w',
                                   'f435w', 'f450w', 'f555w', 'f606w','f814w',
                                   'U2', 'B2', 'V2','R2', 'I2','J2','H2','K2',
                                   'u2','g2','r2','i2','z2','f300w2','f336w2',
                                   'f435w2','f450w2','f555w2','f606w2','f814w2',
                                   'Halpha','FUV','NUV'])
        # Merger Related indicates whether the star is
        # a model with a merger in it.
        astroTable['Merger_Related'] = [False] * len(astroTable)
        # Creates container directory for stellar models grouped by metallicity
        # and whether they are single or binary related
        # the latter included mergers, QHE Models
        if not os.path.isdir(bevo_dir + "iso/" + metallicity + '/' +
                                         "sin" +
                                         "/"):
                        os.makedirs(bevo_dir +
                                    'iso/' + metallicity +
                                    "/" + SinOrBin + "/")
        if not os.path.isdir(bevo_dir + "iso/" + metallicity + '/' +
                                         "bin" +
                                         "/"):
                        os.makedirs(bevo_dir +
                                    'iso/' + metallicity +
                                    "/" + "bin" + "/")
          # Single star models that are not the hmg models go into the <Metallicity>/sin
          # directory.
            if (x[0:10] == 'NEWSINMODS' and 
                not re.matches(r"NEWSINMODS/z[0-9]+_hmg/.*")):
                astroTable.write(bevo_dir + 'iso/' + met + '/' +
                                     "sin" + "/" + "FitsModels" +
                                     x[16:] + ".fits", format='fits',
                                     overwrite=True)
                        print(x[16:] + ".fits")
            else:
                # Models that count as secondary star with compact remnants
                # go into the <metallicity>/bin subdirectory of iso.
                if (x[11:21] == 'NEWSECMODS'):
                            astroTable.write(bevo_dir + 'iso/' + met + '/' +
                                             "bin" + "/" + "FitsModels" +
                                             x[29:] + "sec.fits", format='fits',
                                             overwrite = True)
                        print(bevo_dir + 'iso/' + metallicity + '/' +
                              SinOrBin + "/" + "FitsModels" +
                              x[29:] +"sec.fits")
                else:
                    # Check for merger related (models_type=0) models
                    # and mark them.
                    if (find_mergers(astroTable)):
                        astroTable['Merger_Related'] = [True] * len(astroTable)
                    # Models that count as secondary star with compact remnants
                    # go into the <metallicity>/bin subdirectory of iso.
                    astroTable.write(bevo_dir + 'iso/' + metallicity + '/' +
                                     "bin" + "/" + "FitsModels" +
                                     x[27:] + "bin.fits", format='fits',
                                     overwrite = True)
                    print(x[27:] + ".fits")
                
                except FileNotFoundError:
                    print("File not Readable/Found")
                    notFound.append(x)
                    #The following is commented out since I have yet to delete this version yet.
"""        still_to_do=set()
        
            
            if (x=="zem5" or x=="zem4" or x=="z001" or x=="z002" or x=="z003" or x=="z004"):
                c=set(c).union(set(glob.glob(hoki.MODELS_PATH+"/NEWSINMODS/"+x+"hmg"+"/*")))
            still_to_do=still_to_do.union(a, b, c)
        still_to_do=still_to_do.difference(included)
        lenModPath=len(HOKI.MODELS_PATH)
        for x in still_to_do:
            try:
                f = pd.read_csv(x, sep="\s+",
                                header=None)
                f.apply(pd.to_numeric, errors='coerce')
                f.transpose()
                f.drop(columns=[f.columns[x] for x in colsToElim], inplace=True)
                included.add(hoki.MODELS_PATH+x)
                # Source: https://stackoverflow.com/questions/483666/
                # reverse-invert-a-dictionary-mapping
                invmap = {u: v for v, u in hoki.dummy_dict.items()}
                f.rename(columns=invmap, inplace=True)
                f["model_type"]=[typ]*len(f.index)
                astroTable = Table.from_pandas(f)
                astroTable.remove_columns(['?','modelimf','mixedimf', 'V-I', 'U',
                                           'B', 'V', 'R', 'I', 'J', 'H', 'K','u',
                                           'g', 'r', 'i', 'z', 'f300w', 'f336w',
                                           'f435w', 'f450w', 'f555w', 'f606w','f814w',
                                           'U2', 'B2', 'V2','R2', 'I2','J2','H2','K2',
                                           'u2','g2','r2','i2','z2','f300w2','f336w2',
                                           'f435w2','f450w2','f555w2','f606w2','f814w2',
                                           'Halpha','FUV','NUV'])
                if not os.path.isdir(bevo_dir + "iso/" + metallicity + '/' +
                                     SinOrBin + "/"):
                    os.makedirs(bevo_dir + 'iso/' + metallicity +
                                "/" + SinOrBin + "/")
                if (x[lenModPath+0:lenModPath+10]=='NEWSINMODS'):
                    nameMatcher=re.match(re".*/NEWSINMODS/[a-z0-9]+hmg/.*")
                    if (bool(nameMatcher)):
                        astroTable["model_type"]=[4]*len(astroTable)
                    else:
                        astroTable["model_type"]=[-12]*len(astroTable)
                    astroTable.write(bevo_dir + 'iso/' + met + '/' +
                                     SinOrBin + "/" + "FitsModels" +
                                     x[lenModPath+16:] + ".fits", format='fits',
                                     overwrite=True)
                    print(x[lenModPath+16:]+".fits")
                else:
                    print(x[lenModPath+11:lenModPath+21])
                    
                    if (x[lenModPath+11:lenModPath+21]=='NEWSECMODS'):
                        astroTable.write(bevo_dir + 'iso/' + met + '/' +
                                         SinOrBin + "/" + "FitsModels" +
                                         x[lenModPath+29:] + "sec.fits", format='fits',
                                         overwrite=True)
                        print(bevo_dir + 'iso/' + metallicity + '/' +
                              SinOrBin + "/" + "FitsModels" +
                              x[lenModPath+29:] +"sec.fits")
                        # Added designation for secondary models that can be
                        # either regular or single star
                        astroTable["model_type"]=[-15]*len(astroTable)
                    
                    else:
                        astroTable.write(bevo_dir + 'iso/' + metallicity + '/' + SinOrBin +
                                         "/" + "FitsModels" +
                                         x[lenModPath+27:] +"bin.fits", format='fits',
                                         overwrite=True)
                        # Added designation for secondary models that can be
                        # either regular or single star
                        astroTable["model_type"]=[-10]*len(f.index)
                        print(x[lenModPath+27:] + ".fits")
                        """
                
                except FileNotFoundError:
                    print("File not Readable/Found")
                    notFound.append(x)

    """ The follwing method takes in the reformatted stellar model files listed in
    the corresponding model-metallicity's cluster input files originally meant
    for color magnitude diagram makers. This is suppoesed to create the isochrone
    files out of which BPASS creates its isochrone objects. For each stellar model
    file, we try to find the row in the stellar model representing a model with the
    age closest to the inputted age and is within the specified margin error for
    log(age in years). Then that row is concattenated with other similar rows from
    different stellar model files to create the isochrone table. The age closest 
    method may get rid of some data but is meant to improve accuracy of the age of the
    stellar model and improve computational efficiency to some extent."""        
def extractor(SinOrBin, model, age, metallicity, input_dir, bpass_evo_dir,
              margin):
    entries = glob.glob(input_dir + "/*")
    print(len(entries))
    stopLen = 4
    bigOne = None;
    initlMass = np.nan
    initlMass2 = np.nan
    for x in entries:
        if x[0] != '.' and not os.path.isdir(x):
            # To Obtain initial mass I will use Regex
            # to extract the initial mass of a model.
            try:
                org = Table.read(x, format = 'fits')
                org = org.to_pandas()
                starIndex = "Not Applicable"
                if (x[-8:-5] == 'bin' and org['model_type'][0] == 0):
                    starIndex = findLastRelMax(org)
                f = org[np.where(np.abs(np.log10(org['age'])-age) <= margin)[0]]
                if (len(f) != 0):
                    f=f[np.where(np.abs(f['age'] - 10 ** age) == 
                                 np.min(np.abs(f['age'] - 10 ** age)))[0]]
                if len(f) != 0:
                    index = f.index
                    indexlen = len(index)
                    if initial:
                        initial = False
                        bigOne = f
                        # mass and mass2 stand for the initial masses
                        # of the primary and secondary respectively.
                        bigOne['mass'] = [initlMass] * indexlen
                        bigOne['mass2'] = [initlMass2] * indexlen
                        bigOne['single'] = [True] * indexlen
                        bigOne['mergered?'] = [False] * indexlen
                        if (x[-8:-5] == 'bin' or x[-8:-5] == 'sec'):
                            bigOne['single'] = [False] * indexlen
                            # I only need to check one row since
                            # whether a model is a merger model
                            # is the same for all rows of the model
                            if bigOne['Merger_Related'][0]:
                                # Find the row in which the model merges.
                                # Note that I pass in the original grid.
                                merge_pt = find_mergers(org)
                                bigOne['mergered?'] = [x>merge for x in bigOne.index]
                        print(bigOne.shape)
                    else:
                        f['mass'] = [initlMass] * indexlen
                        f['mass2'] = [initlMass2] * indexlen
                        f['single'] = [True] * indexlen
                        f['Merger Related']=[False] * indexlen
                        if (x[-8:-5] == 'bin'):
                            f['single'] = [False] * indexlen
                            if bigOne['Merger_Related'][0]:
                                merge_pt = find_mergers(org)
                                bigOne['mergered?'] = [x > merge for x
                                                    in bigOne.index]
                        bigOne = pd.concat([f, bigOne], axis = 0)
                        print(bigOne.shape)
            except FileNotFoundError:
                print('File Not Found/Readable')
    
    try:
        if not isinstance(bigOne, type(None)) and not (bigOne['age'].empty):
            bigOne = bigOne.apply(pd.to_numeric, errors = 'coerce')
            reduced = Table.from_pandas(bigOne)
            if not os.path.isdir(bpass_evo_dir +
                                 'iso/' + model + '/' + metallicity + '/'):
                os.makedirs(bpass_evo_dir +
                            'iso/'+ model + '/' + metallicity + '/')
            print(reduced.columns)
            reduced.write(bpass_evo_dir + 'iso/' + model+"/" + metallicity +
                          '/' + 'iso' + str(age) + SinOrBin +
                           '.fits', format = 'fits', overwrite = True)
            print(bigOne.shape)
    except IndexError:
        print ('It looks like there are no stars in out input file' +
               ' that have the specified age...')