#Project
import numpy as np
import pylab as plt
import math

def make_star(x_axis,y_axis,x_centre,y_centre):
    #Makes a stellar model using the
    #axis lengths and centre points
    size=x_axis*y_axis           #Initialise size
    R=(x_axis/2)-20              #Radius of star with buffer of 20
    #Create lists for x- and y- axes
    x=np.arange(size*1.) % x_axis
    y=np.arange(size*1.) % y_axis    
    #Reshape to form arrays
    x = np.reshape(x, (x_axis,y_axis))
    y=np.rot90(np.reshape(y, (x_axis,y_axis)))
    #Subtract centre points
    x=x-x_centre
    y=y-y_centre    
    #Square array
    x=np.square(x)
    y=np.square(y)    
    r1=np.add(x,y)
    #Square root arrays to obtain circles
    r=np.sqrt(r1)

    #Make circle of stellar radius
    for i in range (0,(x_axis-1)):
        for j in range(0,(y_axis-1)):
            #Outside the stellar radius will be zeroed     
            if (r[i][j]>R):
                r[i][j]=0.0
    #Plot stellar model    
    plt.imshow(r, cmap=plt.cm.binary)

def limb_darkening(a,b,gamma):
    #Calculates quadratic limb darkening coefficient
    mu=math.cos(gamma)
    LD=1-(a*(1-mu))-(b*((1-mu)**2))
    print (LD)
    
make_star(1024,1024,512,512)
limb_darkening(0.3,0.1,0)
