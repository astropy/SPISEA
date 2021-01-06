# Makes the Isochrone Files
from spisea import IsochroneMakerReformattedVersion2
from spisea.IsochroneMakerReformattedVersion2 import extractorO, reformatterO
from spisea import IsochroneMakerReformattedVersionNew
from spisea.IsochroneMakerReformattedVersionNew import reformatter, extractor
import time
import math
# Note: O exists in order to confuse between the extractors and reformatters from the two functions
reformatterO("/g/lu/scratch/ryotainagaki/BPASS_tester_old_style/" ,["z020"])
for x in range(30):
    t1 = time.time()
    extractorO(round(7.0 + 0.1 * x,1),"z020","/g/lu/scratch/ryotainagaki/" +
               "BPASS_tester_old_style/",
               "/g/lu/scratch/ryotainagaki/BPASS_iso_filesTimed/", 0.05)
    t2 = time.time()
    print(t2 - t1)
reformatter("/g/lu/scratch/ryotainagaki/BPASS_tester_newReformatTest/", ["z020"])
print("Begin New Method")
for x in range(30):
    t1 = time.time()
    extractor(round(7.0 + 0.1 * x,1),"z020","/g/lu/scratch/ryotainagaki/" +
              "BPASS_tester_newReformatTest/",
              "/g/lu/scratch/ryotainagaki/BPASS_iso_filesTimedIsolated/", 0.05)
    t2 = time.time()
    print(t2 - t1)