#make 1D vertical array
import numpy as np

size = 1024
planet=np.arange(0,size*1.)
centre = size/2
planet[centre]=1

while i < size:
    planet