# Makes the Isochrone Files
from spisea import IsochroneMakerReformattedVersionNew
from spisea.IsochroneMakerReformattedVersionNew import reformatter, extractor
import time
import math

for met in ["zem5", "zem4", "zem3", "z001", "z002", "z003", "z004", "z005", "z008", "z010", "z014", "z020", "z030", "z040"]:
    for x in range(51):
        t1=time.time()
        # TODO: replace "/g/lu/scratch/ryotainagaki/BPASS_tester_newReformatTest/" with absolute path to your directory of
        # reformatted fitsmodels
        # that were created by the reformatter function
        # TODO: replace "/g/lu/scratch/ryotainagaki/BPASS_iso_filesTimedIsolated/" with absolute path to your intended directory
        # of isochrone files
        extractor(round(6.0+0.1*x,1),met,"/g/lu/scratch/ryotainagaki/BPASS_tester_newReformatTest/","/g/lu/scratch/ryotainagaki/BPASS_iso_filesTimedIsolated/",0.05)
        t2=time.time()
        print(t2-t1)