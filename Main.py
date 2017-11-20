from Lib.ParticleMove import ParticleMove
from Lib.RepresentData import *
import numpy as np
from Lib.Functions import *
import time
from random import random
import datetime

d = ReadXml()

E_listWOS = SetupElementsWOS(d)

EWOS = ResultingString(E_listWOS)

t1 = time.time()
print(datetime.datetime.now().time())
x = 1.25  # ook bv x = 3*y
y = np.arange(-3, 3.05, 0.1).tolist()
z = np.arange(0, 2.05, 0.1).tolist()


Vgrid = WalkOnSpheres_potential_slice(EWOS,x,y,z)
t2 = time.time()
print(t2-t1)

print(Vgrid.shape)




#https://www.howtoforge.com/tutorial/install-git-and-github-on-ubuntu-14.04/
