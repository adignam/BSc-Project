import numpy as np
import glob as glob
import os
from astropy.io import fits

#Specify path and file list
path = r'\\ad\data\STUHOME\13\43\dignama3\Python Scripts\ccf'
filelist = glob.glob(os.path.join(path, '*.*') )
JD = np.zeros(len(filelist))

#Use loop to find Julian date for all files in path
i = 0
for infile in filelist:
    
    hdulist = fits.open(infile)
    MJD = hdulist[0].header['mjd-obs']
    v1=hdulist[0].header['CRVAL1']
    v2=hdulist[0].header['CDELT1']
    v3=hdulist[0].header['NAXIS1']
    v=(v1+v2)*np.arange(0,v3)
    JD[i] = (MJD + 2400000.5)
    MJDEQ = str(JD[i])
    Entry = str(i) 
    print("Current file number " + Entry + " is: " + infile + " and its Julian date is " + MJDEQ)
    hdulist.close()
    i = i + 1
    
#Display JD array
print(v)
#Save array into a text file
np.savetxt('JulianDates',JD)

def get_phase():
    path = r'\\ad\data\STUHOME\13\43\dignama3\Python Scripts\ccf'

    filelist = glob.glob(os.path.join(path, '*.*') )
    phase = np.zeros(len(filelist))
    Ins = filelist.copy()
    Info = open('Info.txt','w')
    print("File name                             Phase      Orbit Number   Instrument")
    
    i = 0
    Data = np.loadtxt('JulianDates')
    FinalP=np.zeros(len(filelist))
    for infile in filelist:
            
        hdulist = fits.open(infile)
        Ins[i] = hdulist[0].header['HIERARCH ESO INS ID']  
        
        #Calculation for Phase and Orbit Number
        f = filelist[i][36:]
        T = 2454279.436714
        P = 2.21857567
        phase[i] = (Data[i] - T)/P 
        OrbitNo = round(phase[i],0)
        FinalP[i] = phase[i] - OrbitNo
        
        #print(("(%s) %3s %10f %7i %15s" % (i,f,FinalP,OrbitNo,Ins[i])),file = Info)
        #print("(%s) %3s %10f %7i %15s" % (i,f,FinalP,OrbitNo,Ins[i]))
        hdulist.close()
        i += 1
    return FinalP
get_phase()
