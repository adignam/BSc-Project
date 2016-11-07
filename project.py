import numpy as np
import pylab
import math

size = 1025
r=512
x=0
y=0

star=np.zeros((size,size))

for x in range (0,size):
    for y in range (0,size):
        dist = math.sqrt(abs(pow((x-r),2)+pow((y-r),2)))
        if dist < r :
            star[x,y]=1
        y+=1
    x+=1

print (star)

pylab.imshow(star, origin='lower')