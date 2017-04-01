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

def gauss(a,sigma,v):
    vrot=(np.arange(1024)-512)/512.*2.5
    #v=np.arange(-12,12,0.02345)   
    x=np.shape(vrot)
    y=np.shape(v)
    x=np.rint(x)
    y=np.rint(y)
    gauss=np.zeros((1024,161))
    i=0
    for vx in vrot:
        gauss[i,:]=1-a*np.exp((-(v-vx)**2)/(2*sigma**2))
        i=i+1
    return gauss
    
def spectrum(array,g,v):  
    #Plots the spectrum of the star
    
    flux=np.zeros(np.shape(v))
    fstar=np.sum(array,0) 
    for i in range(0,161):       
        flux=flux+fstar[i]*g[i,:]
        i+=1
    flux=flux/np.median(flux[((v<-12)|(v>12))])        
    return flux   
    
    
def planet_motion(aRstar,incl,Rstar,phase,phase_bin,centre,v,amp,line_width):
    #phase=np.arange(-phase,phase,phase_bin)
    x=aRstar*np.sin(2*np.pi*phase)*Rstar+centre   
    y=aRstar*np.cos(incl)*np.cos(2*np.pi*phase)*Rstar+centre
    Rplanet = Rstar/10.
    g=gauss(amp,line_width,v)
    
    star,mask1=make_circle(1025,1025,Rstar)
    star=limb_darkening(0.3,0.3,star/Rstar,mask1)
    star=np.nan_to_num(star)  
    print(np.shape(v))
    flux=np.zeros( (111,161) )
    for i in range(0 ,len(phase)):
        planet,mask=make_circle(1025,1025,Rplanet,x[i],y[i])
        system=star*(1.-mask)
        flux[i,:]=spectrum(system,g,v)
        #plt.imshow(system, cmap=plt.cm.afmhot)
        #plt.title(i)        
        #plt.plot(v,flux[i,:])
        #plt.pause(0.5)
        #plt.close()
    #plt.imshow(flux, cmap=plt.cm.afmhot)    
    return flux

def get_phase():
    path = r'/home/member/a/aishling/Project/ccf'
    filelist = glob.glob(os.path.join(path, '*A.*') )
    phase0 = np.zeros(len(filelist))
    phase=np.zeros(len(filelist))
    data=np.zeros((len(filelist),161))
    i = 0
    T = 2454279.436714
    P = 2.21857567
    
    for infile in filelist:            
        hdulist = fits.open(infile)
        MJD = hdulist[0].header['mjd-obs']
        Date = (MJD + 2400000.5)        
        phase0[i] = (Date - T)/P 
        OrbitNo = round(phase0[i],0)
        phase[i] = phase0[i] - OrbitNo

        v1=hdulist[0].header['CRVAL1']
        v2=hdulist[0].header['CDELT1']
        v3=hdulist[0].header['NAXIS1']
        v_data=v1+v2*np.arange(v3)+2.55

        im=hdulist[0].data
        spec=im[-1,:]
        spec=spec/np.median(spec[((v_data<-12)|(v_data>12))])
        data[i,:]=spec
        hdulist.close()
        i += 1
    plt.imshow(data)
    return phase, data, v_data

def mcmc(parameters,stepsize,phase,data,v):
    constant=[0,1,2,3,4,5]
    constant[0]=1.           #8.84aRstar
    constant[1]=85.71*np.pi/180.   #incl
    constant[2]=500               #Rstar
    constant[3]=512              #centre
    constant[4]=0.816                   #u1
    constant[5]=0                       #u2

    nstep=200
    npars=np.size(parameters)                          
    error=0.00015
    u1=0.3
    u2=0.3
    #phase=0.1
    chain=np.zeros((npars,nstep))
    chi_squared=np.zeros(nstep)
    chain[:,0]=parameters

    model=planet_motion(constant[0],constant[1],constant[2],phase,0.001,constant[3],v,parameters[2],parameters[3])
    flux=data

    chisq=np.sum((flux-model)**2/error**2 )
    chi_squared[0]=chisq
    
    for i in range(1,nstep):
        ##Make a random jump in parameter space
        new_parameters=chain[:,i-1]+np.random.normal(0.,1.,npars)*stepsize
        model=planet_motion(constant[0],constant[1],constant[2],phase,0.001,constant[3],new_parameters[4],new_parameters[2],new_parameters[3])
        nchisq=np.sum((flux-model)**2/error**2 )
        dchi=nchisq-chisq
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
    print(new_parameters)
    return chain,chi_squared

Rp_Rstar=0.151
lamda=-0.85*np.pi/180. ##radians
amplitude=0.9
line_width=2.7
vsini=3.5

e_Rp_Rstar=0.00012
e_lamda=0.28
e_amplitude=0.02
e_line_width=0.02
e_vsini=0.2
stepsize=[e_Rp_Rstar,e_lamda,e_amplitude,e_line_width,e_vsini]
phase,data,v=get_phase()
parameters=[Rp_Rstar,lamda,amplitude,line_width,vsini]
mcmc(parameters,stepsize,phase,data,v)
