import os
import glob
from astropy.table import Table

def convert():

    """
    Converts evolution model .dat files to fits tables. Saves new files with the same
    base name, but .fits extension, in the original folder.

    At runtime, will ask for path to the .dat files, plus the file prefix. This is
    generally 'iso' or something similar, but varies. Enter in both without quotes.
    """
    print "Enter the isochrone path (no quotes): " ; inpath = raw_input()
    inpath=inpath.strip()
    print "Enter the isochrone file prefix (no quotes): "; fileprefix = raw_input()
    fileprefix=fileprefix.strip()

    filelist = glob.glob(inpath+'/'+fileprefix+'*.dat')
    
    for file in filelist:        
        current_data = Table.read(file,format='ascii')
        base = os.path.splitext(file)[0]
        if os.path.isfile(base+'.fits'):
            print "The output exists already. Delete it?"; answer = raw_input()
            if answer in ('y','yes','ok'):
                os.remove(base+'.fits')
            else:
                print "Output already exists. Exiting."
                return
        current_data.write(base+'.fits',format='fits')
