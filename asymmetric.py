#Project
import numpy as np
import pylab as plt

def make_star(x_axis,y_axis):

    size=x_axis*y_axis           #Initialise size
    #radius=(size/2)              #Radius of star
    #add buffer  
    x_centre=x_axis/2
    y_centre=y_axis/2
    #Create list
    x=np.arange(size*1.) % x_axis
    y=np.arange(size*1.) % y_axis
    #print (x)
    #Reshape to form array
    x = np.reshape(x, (x_axis,y_axis))
    #print (x)
    #y array
    y=np.rot90(np.reshape(y, (x_axis,y_axis)))
    #Subtract centre points
    x=x-x_centre
    y=y-y_centre
    #print (x)
    #print (y)
    #Square array
    x=np.square(x)
    #print (x)
    y=np.square(y)
    #Add x and y
    if (x_centre>y_centre):
        i=0
        for i in range (0,y):
            r1=x[i]+y[i]
            i+=1
    arr1 = x[(y_axis-1):int(x_axis-1)]
    print (arr1) 
    #Square root arrays
    #r=np.sqrt(r1)
    #plt.imshow (r)
    
make_star(1024.,1024)