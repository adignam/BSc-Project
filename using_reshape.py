#Project
import numpy as np
import pylab as plt

size=50*50                   #Initialise size
radius=(size/2)              #Radius of star 
i=0
j=0
#Create list
x=np.zeros((size))
print (x)
#Reshape to form array
x = np.reshape(x, (50,50))
print (x)
#Copy array
y=np.rot90(x)
#Populate arrays
for i in range (0,49):
    if (i < radius):
        x[i]=1
    i+=1

for j in range (0,49):
    if (j < radius):
        x[j]=1
    j+=1
    
print (x)
print (y)
#Square array
x=np.square(x*x)
print (x)
y=np.square(y*y)
#Add arrays
r=x+y
print (r)
#Square root arrays
r=np.sqrt(x+y)
plt.imshow (r)
