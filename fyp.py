#Project
import numpy as np
import matplotlib.pyplot as plt
from random import random
import os
import glob
from astropy.io import fits


def make_circle(x_axis,y_axis,R=510,xc=512.,yc=512.):
    size=x_axis*y_axis           #Initialise size
    x_centre=xc#x_axis/2
    y_centre=yc#y_axis/2   
    #Create list
    x=np.arange(size*1.) % x_axis
    y=np.arange(size*1.) % y_axis    
    #Reshape to form array
    x = np.reshape(x, (x_axis,y_axis))
    #y array
    y=np.rot90(np.reshape(y, (x_axis,y_axis)))
    #Subtract centre points
    x=x-x_centre
    y=y-y_centre    
    #Square array
    x=np.square(x)
    y=np.square(y)    
    r1=np.add(x,y)
    #Square root arrays
    r=np.sqrt(r1)
    #Make circle of given stellar radius
    mask=(r<=R)
    r=r*mask
    #Plot stellar model    
    #plt.imshow(mask, cmap=plt.cm.binary)
    #plt.ylim([0,1024])
    return r, mask
    
def limb_darkening(a,b,r,mask):
    #Calculates limb darkening for the star    
    mu=np.sqrt(1-r**2)
    LD=r*0.
    LD[mask==1]=1-(a*(1-mu[mask==1]))-(b*((1-mu[mask==1])**2))
    return LD

def gauss(a,sigma,veq):
    vrot=(np.arange(1024)-512)/512.*veq
    v=np.arange(-12,12,0.02345)  
    
    gauss=np.zeros( (vrot.shape[0],v.shape[0]) )
    i=0
    for vx in vrot:
        gauss[i,:]=1-a*np.exp((-(v-vx)**2)/(2*sigma**2))
        i=i+1
    return gauss,v
    
def spectrum(array,g,v):  
    #Plots the spectrum of the star
    
    flux=np.zeros(np.shape(v))
    fstar=np.sum(array,0) 
    for i in range(0,len(v)):       
        flux=flux+fstar[i]*g[i,:]
        i+=1
    return flux   
    
    
def planet_motion(aRstar,incl,Rstar,phase,phase_bin,centre):
    phase=np.arange(-phase,phase,phase_bin)
    x=aRstar*np.sin(2*np.pi*phase)*Rstar+centre   
    y=aRstar*np.cos(incl)*np.cos(2*np.pi*phase)*Rstar+centre
    Rplanet = Rstar/10.
    g,v=gauss(0.9,2.7,2.)
    
    
    star,mask1=make_circle(1025,1025,Rstar)
    star=limb_darkening(0.3,0.3,star/Rstar,mask1)
    star=np.nan_to_num(star)    
    g,v=gauss(0.9,0.7,2.)
    flux=np.zeros( (len(phase),len(v)) )
    for i in range(0 ,len(phase)):
        planet,mask=make_circle(1025,1025,Rplanet,x[i],y[i])
        system=star*(1.-mask)
        flux[i,:]=spectrum(system,g,v)
        #plt.imshow(system, cmap=plt.cm.afmhot)
        #plt.title(i)        
        #plt.plot(v,flux[i,:])
        #plt.pause(0.5)
        #plt.close()
    plt.imshow(flux, cmap=plt.cm.afmhot)    
    return flux, phase

def get_phase():
    path = r'C:\Users\Aishling\Documents\Python Scripts'
    filelist = glob.glob(os.path.join(path, '*.*') )
    phase = np.zeros(len(filelist))
    Ins = filelist.copy()
    Info = open('Info.txt','w')
    print("File name                             Phase      Orbit Number   Instrument")
    
    i = 0
    Data = np.loadtxt('JulianDate')
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
        
        print(("(%s) %3s %10f %7i %15s" % (i,f,FinalP,OrbitNo,Ins[i])),file = Info)
        print("(%s) %3s %10f %7i %15s" % (i,f,FinalP,OrbitNo,Ins[i]))
        hdulist.close()
        i += 1
    return FinalP
        


def mcmc(parameters,phase,data,v):
    aRstar=1.           #8.84
    incl=90.*np.pi/180.                #85.71
    Rstar=500
    centre=512
    phase=get_phase()
    error=0.3
    stepsize=1
    nstep=1
    npars=np.size(parameters)
    lamda=-0.85           #degrees
    Rplanet=216/Rstar     #normalised by Rstar
    vsini=3.5                             
    period=2.21857567*86400 #seconds
    Mstar=0.806         #solar masses
    
    ##free parameters
    #rprstar
    #veq
    #linewidth (FWHM/2.35)
    #lambda
    #amplitude
    
    ##get v from the fits files
    
    ##fix
    #aRstar
    #incl
    #u1
    #u2
    
    
    flux=planet_motion(phase,aRstar,incl,Rstar,0.1,0.001,centre)
    phase=get_phase()

    #define the array that stores all the steps.
    chain=np.zeros((npars,nstep))
    chi_squared=np.zeros(nstep)
    chain[:,0]=parameters

    ##calculate the model 
    model,phase=planet_motion(aRstar,i,Rstar,0.1,0.001,centre)
    
    ##calculate Chi squared.
    chisq=np.sum( (flux-model)**2/error**2 )
    chi_squared[0]=chisq
    
    ## Now run the MCMC algorithm Nstep times
    for i in range(1,nstep):
        ##Make a random jump in parameter space
        new_parameters=chain[:,i-1]+np.random.normal(0.,1.,npars)*stepsize

        ## make model, calculate current chi squared and compare to previous value
        model=make_circle(1025,1025)
        nchisq=np.sum( (flux-model)**2/error**2 )
        dchi=nchisq-chisq
        print(chi_squared)
        ## Going to lower chi squared: Accept new value
        if dchi <=0.:
            chain[:,i]=new_parameters
            chisq=nchisq
        else:
            rn=np.random.uniform()  ##There is still a chance to go to a higher chi squared...
            p=np.exp(-(dchi)/2.) 
            if rn <= p:
                chain[:,i]=new_parameters   #accepted higher chi squared, move to new position
                chisq=nchisq
            else:
                chain[:,i]=chain[:,i-1]  ## rejected. Stay at current point
        chi_squared[i]=chisq
        
    return chain,chi_squared

#star=make_circle(1024,1024,510)           make star 
#planet,m=make_circle(1024,1024,256)         #make planet
#flux,v=planet_motion(1.,1.,500,0.1,0.001,512) 
#printflux)           # coeff=8.84,i = 1.14959 for hd189733b
Rp_Rstar=0.12
lamda=-0.85  *np.pi/180. ##radians
amplitude=0.9
line_width=2.7
vsini=2.

e_Rplanet=0.027
e_lamda=+0.28

parameters=[Rplanet,lamda,amplitude,line_width,vsini]
mcmc(parameters)
#plt.imshow(planet, cmap=plt.cm.binary)
#scipy.misc.imsave('project_image.jpg', -star)
