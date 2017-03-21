#Project
import numpy as np
import matplotlib.pyplot as plt
from IPython import display
import math
import scipy.misc
import time

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
    for vx in vrot:
        gauss=1-a*np.exp((-(v-vx)**2)/(2*sigma**2))
    return gauss,v
    
def spectrum(r,v,gauss):  
    #Plots the spectrum of the star
    flux=np.zeros(np.shape(v))
    fstar=-np.sum(r,0) 
    for i in range(0,len(v)):       
        flux=fstar[i]*gauss
        i+=1
    return flux    
    
def planet_motion(coeff,i,Rstar,phase,phase_bin,centre):
    phase=np.arange(-phase,phase,phase_bin)
    x=coeff*np.sin(2*np.pi*phase)*Rstar+centre   
    y=coeff*np.cos(i)*np.cos(2*np.pi*phase)*Rstar+centre
    Rplanet = Rstar/2
    
    star,mask1=make_circle(1025,1025,Rstar)
    star=star-limb_darkening(0.3,0.3,star,~mask1)
    star=np.nan_to_num(star)

    g,v=gauss(0.9,0.7,2.)

    for i in range(0,200):
        if (x[i]<0):
            planet,mask=make_circle(1025,1025,Rplanet,-x[i],y[i])
        else:
            planet,mask=make_circle(1025,1025,Rplanet,x[i],y[i])
        planet=limb_darkening(0,0,planet,mask)       
    plt.figure()
    plt.title(i)
    mask=1-mask
    system=star+planet
    #plt.imshow(system, cmap=plt.cm.binary)
        #plt.close()
    flux=spectrum(system,v,g)
    plt.plot(v,-flux)
 
#star=make_circle(1024,1024,510)           make star 
#planet,m=make_circle(1024,1024,256)         #make planet
planet_motion(1,1.,500,0.1,0.001,512)            # coeff=8.84,i = 1.14959 for hd189733b

#plt.imshow(planet, cmap=plt.cm.binary)
#scipy.misc.imsave('project_image.jpg', -star)

==============================================
#Project
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.misc
from random import random

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
    for vx in vrot:
        gauss=1-a*np.exp((-(v-vx)**2)/(2*sigma**2))
    return gauss,v
    
def spectrum(r,v,gauss):  
    #Plots the spectrum of the star
    flux=np.zeros(np.shape(v))
    fstar=-np.sum(r,0) 
    for i in range(0,len(v)):       
        flux=fstar[i]*gauss
        i+=1
    return flux    
    
def planet_motion(coeff,i,Rstar,phase,phase_bin,centre):
    phase=np.arange(-phase,phase,phase_bin)
    x=coeff*np.sin(2*np.pi*phase)*Rstar+centre   
    y=coeff*np.cos(i)*np.cos(2*np.pi*phase)*Rstar+centre
    Rplanet = Rstar/2
    
    star,mask1=make_circle(1025,1025,Rstar)
    star=star-limb_darkening(0.3,0.3,star,~mask1)
    star=np.nan_to_num(star)

    #g,v=gauss(0.9,0.7,2.)

    for i in range(0,200):
        if (x[i]<0):
            planet,mask=make_circle(1025,1025,Rplanet,-x[i],y[i])
        else:
            planet,mask=make_circle(1025,1025,Rplanet,x[i],y[i])
        planet=-limb_darkening(0,0,planet,mask)   
        plt.figure()
        plt.title(i)
        system=star+planet
        plt.imshow(system, cmap=plt.cm.binary)        
        #print(x[i],y[i])
    #flux=spectrum(star,v,g)
        
    
def chi_squared(data,model,sigma,v):
    sigma=1.        
    v1=[]
    v2=[]
    chi1=sum((data-model)**2)/sigma**2
    chi2=sum((data-model)**2)/sigma**2
    
    for i in v:
        v1.append(v1)
        v2.append(v2)

        chi1.append(chi1)
        chi2.append(chi2)
        
        v2[i]=v1[i]+random(1)        
    
    ratio=np.exp(((chi1-chi2)/2)*sigma**2)
    
    if random() < ratio:
        v1=v2
        chi1=chi2
    
        
#star=make_circle(1024,1024,510)           #make star 
#planet,m=make_circle(1024,1024,256)         #make planet
planet_motion(1,1.,512,0.1,0.001,512)            # coeff=8.84,i = 1.14959 for hd189733b

#plt.imshow(planet, cmap=plt.cm.binary)
#scipy.misc.imsave('project_image.jpg', -star)
