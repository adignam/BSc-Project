#Project
import numpy as np
import pylab as plt
import math

def make_star(x_axis,y_axis):
    size=x_axis*y_axis           #Initialise size
    R=(x_axis/2)-20              #Radius of star with buffer of 20
    x_centre=x_axis/2
    y_centre=y_axis/2
    
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
    #plt.imshow(r, cmap=plt.cm.binary)
    #plt.ylim([0,1024])
    return r
    
def limb_darkening(a,b,gamma):
    mu=math.cos(gamma)
    LD=1-(a*(1-mu))-(b*((1-mu)**2))
    return LD
    
def spectrum(r, A, sigma):  
    #Equatorial Velocity
    veq=2.
    vrot=(np.arange(1025)-512)/512.*veq
    v=np.arange(-12,12,0.02)  #-12....12
    flux=np.zeros(1025)

    #Gaussian
    fstar=np.sum(r,0)
    i=0
    for vx in vrot :
        #Gaussian calculation
        Gauss=1-A*np.exp((-(v-vx)**2)/(2*sigma**2))
        print("Gauss=",Gauss)
        #Calculate flux
        flux[i,:]=fstar[i]*Gauss
        i+=1
    plt.plot(v,flux.sum(1))
    
    
    
r=make_star(1024,1024)
r=r-limb_darkening(0.3,0.1,1)
spectrum(r, 0.9, 2)
#plt.imshow(r, cmap=plt.cm.binary)
