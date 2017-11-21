from Lib.ParticleMove import ParticleMove
from Lib.RepresentData import *
import numpy as np
from Lib.Functions import *
import time
from random import random
import datetime

d = ReadXml()


#P_var reset niet

"""print(d)

P = [5, 6, 1]

dist1 = dist_cylinder(P, [1.25, 1.25, 0.3], [1.25, 1.25, 1.3], 0.5)
dist2 = dist_circdisc_4holes(P, [0, 0, 0], [1.25, 1.25, 0], [-1.25, -1.25, 0], [-1.25, 1.25, 0], [1.25, -1.25, 0], 3000, 0.8, 0.8, 0.8, 0.8, [0, 0, 1])

print(dist1)
print(dist2)"""

E_listWOS = SetupElementsWOS(d)

EWOS = ResultingString(E_listWOS)


print(np.zeros(shape=(5,5,5)))

space = 10 ** (-2)
dim_y = np.arange(-10, 10.01, space).tolist()
dim_y = np.arange(-10, 10.01, space).tolist()
dim_z = np.arange(0, 5.01, space).tolist()













"""
#dimensions of grid where all the minimum distances will come in:



t1 = time.time()
print(datetime.datetime.now().time())
x = 0  # ook bv x = 3*y
y = np.arange(-3, 3.05, 0.4).tolist()
z = np.arange(0, 2.05, 0.4).tolist()

y = [1.4]
z = [1.5]

Vgrid = WalkOnSpheres_potential_slice(EWOS,x,y,z)
t2 = time.time()
print(t2-t1)"""





#https://www.howtoforge.com/tutorial/install-git-and-github-on-ubuntu-14.04/
#cd /home/robin/Documents/code
#git add *
#git commit -m "ee"
#git push origin master
