#Find Modified Julian date for each file in a folder
from astropy.io import fits
import os
import glob
import numpy as np

#Specify path and file list
path=r'C:\Users\Aishling\Documents\Python Scripts\ccf'
filelist = glob.glob(os.path.join(path, '*.*') )
JD = np.zeros(len(filelist))

#Use loop to find Julian date for all files in path
i = 0
for infile in filelist:
    
    hdulist = fits.open(infile)
    MJD = hdulist[0].header['mjd-obs']
    JD[i] = (MJD + 2400000.5)
    MJDEQ = str(JD[i])
    Entry = str(i) 
    print("Current file number " + Entry + " is: " + infile + " and its Julian date is " + MJDEQ)
    hdulist.close()
    i = i + 1
    
#Display JD array
print(JD)
#Save array into a text file
np.savetxt('JulianDates',JD)