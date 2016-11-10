#Project
import numpy as np

size=50*50                   #Initialise size 1048576
radius=(size/2)
i=0
#Create list
x=np.zeros((size))
print (x)

#Reshape to form array
x = np.reshape(x, (50,50))
print (x)

#Copy array
y=np.rot90(x)

for i in range (0,size):
    if (x[i] < radius):
        x[i]=1
    i+=1

print (x)
    

'''
#Square array
x*x
y*y

#Add arrays
x+y
#Square root arrays
sqrt(x+y)
'''