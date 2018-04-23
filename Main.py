# main function from which everything in the Lib will be called and the program will be executed

from Lib.Particle import *
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
                                                            # (some electrical fields will be calculated by the WOS method (see functions_WOS) in the particle class itself)

if d["MagneticFieldPlot"] == "yes":                         # plot the magnetic field if enabled in the xml-file
    Plotfield(B_analytic, "magnetic", d)

if d["ElectricFieldPlot"] == "yes":                         # plot the electric field if enabled in the xml-file
    Plotfield(E_analytic, "electric", d)

trajectory = particle.ParticleMove(B_analytic, E_analytic, electrodes, electrodes_WOS, d)               # calculate trajectory
if d["WriteDataToFile"] == "yes":
    particle.WriteToFile(E_analytic, B_analytic, trajectory, d)                                         # write data concerning the trajectory of the particle to a file if enabled in the xml-file
if d["TrajectoryPlot"] == "yes":
    particle.PlotTrajectory(trajectory)                                                                 # plot the trajectory if enabled in the xml-file









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