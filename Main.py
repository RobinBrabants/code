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



electrodes, electrodes_WOS, particle, d = ReadXml()         # extract all the data from the xml file

B_analytic, E_analytic = ResultingField(electrodes)         # already calculates the analytical fields for those objects for which it is possible
                                                            # (some electrical fields will be calculated by the WOS method in the particle class itself)


data = ParticleMove(B, E, "EWOS", d)

if d["WriteDataToFile"] == "yes":
    WriteFile(B, data, d)

if d["MagneticFieldPlot"] == "yes":
    Plotfield(B_analytic, d)

if d["ElectricFieldPlot"] == "yes":
    Plotfield(E_analytic, d)

if d["TrajectoryPlot"] == "yes":
    PlotTrajectory(data)











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