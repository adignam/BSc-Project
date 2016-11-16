#Project
import numpy as np
import pylab as plt
import math

def make_star(x_axis,y_axis):
    size=x_axis*y_axis           #Initialise size
    #add buffer  
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
    plt.imshow (r)
    
def limb_darkening(a,b,gamma):
    mu=math.cos(gamma)
    LD=1-(a*(1-mu))-(b*((1-mu)**2))
    print (LD)
    
make_star(1024.,1024.)
limb_darkening(0.3,0.1,57)