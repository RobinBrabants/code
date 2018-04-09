# main function from which everything in the Lib will be called

from Lib.ParticleMove import ParticleMove
from Lib.RepresentData import *
import numpy as np
from Lib.Functions import *
import time
from random import random
import datetime
from Lib.Elements import *
import xml.etree.ElementTree as ET
import inspect
import sys





electrodes_WOS = []




B_list = SetupElements(d)

B = ResultingField(B_list)

E_list = SetupElements2(d)

E = ResultingField(E_list)



#WalkOnSpheres(EWOS)


data = ParticleMove(B, E, "EWOS", d)

if d["WriteDataToFile"] == "yes":
    WriteFile(B, data, d)

if d["MagneticFieldPlot"] == "yes":
    Plotfield(B, d)

if d["ElectricFieldPlot"] == "yes":
    Plotfield(E, d)

if d["TrajectoryPlot"] == "yes":
    PlotTrajectory(data)







"""
WOS:


d = ReadXml()

print(datetime.datetime.now().time())
x = 5  # ook bv x = 3*y
y = np.arange(-4, 4, 2).tolist()
z = np.arange(-4, 4, 2).tolist()


t1 = time.time()
Vgrid = WalkOnSpheres_potential_slice(d,x,y,z)

t2 = time.time()

print(t2-t1)
"""





"""
https://www.howtoforge.com/tutorial/install-git-and-github-on-ubuntu-14.04/
cd /home/robin/Documents/code
git add *
git commit -m "ee"
git push origin master
"""
#username: RobinBrabants

#http://chris35wills.github.io/conda_python_version/

#https://pythonhosted.org/PyInstaller/usage.html


#packages:
"""
python -m pip install scipy
python -m pip install sympy
python -m pip install matplotlib
apt-get install python3-tk
"""