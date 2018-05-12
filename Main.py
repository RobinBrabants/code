# main function from which everything in the Lib will be called and the program will be executed


from Lib.Particle import *
from Lib.RepresentData import *
from Lib.Functions import *
import sys

import time
import multiprocessing as mp



""""# TEST WOS:
from Lib.Functions_WOS import *
from sympy.vector import CoordSys3D
import sympy as sy

Point = [3, 2, 10]

L = CoordSys3D('L')
x, y, z = sy.symbols('x y z')

electrodes, electrodes_WOS, particles, d = ReadXml()        # extract all the data from the xml file
print(electrodes)
print(electrodes_WOS)

B_analytic, E_analytic = ResultingField(electrodes)
E, B = GetFields(Point, Vector.zero, B_analytic, E_analytic, [], d)

print("\r\nanalytical: " + str(E))
E_WOS = ElectricalField_WOS(electrodes_WOS, Point, d)
print("WOS: " + str(E_WOS))

E = E.components
E_WOS = E_WOS.components
print("WOS/analytical:")
print(E_WOS[L.i]/E[L.i], E_WOS[L.j]/E[L.j], E_WOS[L.k]/E[L.k])"""



electrodes, electrodes_WOS, particles, d = ReadXml()        # extract all the data from the xml file


for particle in particles:
    # check whether the starting position of the particle lies out of the predetermined box
    if particle.Position[0] < d["xmin"] or particle.Position[0] > d["xmax"] or particle.Position[1] < d["ymin"] or particle.Position[1] > d["ymax"] or particle.Position[2] < d["zmin"] or particle.Position[2] > d["zmax"]:
        sys.exit("ERROR: The Starting position of %s (%s) lies out of the predetermined box, check your xml-file\r\n" % (particle.name, particle.Type))

    # check if the starting position of the particle is inside an object
    for electrode in electrodes + electrodes_WOS:
        if electrode.IsPointInObject(particle.Position, d["interval"]) == 2:
            sys.exit("ERROR: The Starting position of %s (%s) lies inside %s, check your xml-file\n" % (particle.name, particle.Type, electrode.name))


print("calculating analytical fields...")
B_analytic, E_analytic = ResultingField(electrodes)         # already calculates the analytical fields for those objects for which it is possible







for particle in particles:
    B_analytic, E_analytic = ResultingField(electrodes)
    str(B_analytic)
    #print(B_analytic)
    trajectory = particle.ParticleMove(B_analytic, E_analytic, electrodes, electrodes_WOS, d)               # calculate trajectory

    if d["WriteDataToFile"] == "yes":
        particle.WriteToFile(E_analytic, B_analytic, d)                                         # write data concerning the trajectory of the particle to a file if enabled in the xml-file
    if d["TrajectoryPlot"] == "yes":
        particle.PlotTrajectory()                                   # plot the trajectory if enabled in the xml-file





"""# multiprocess test:
def multiprocess(particle):
    B_analytic, E_analytic = ResultingField(electrodes)
    trajectory = particle.ParticleMove(B_analytic, E_analytic, electrodes, electrodes_WOS, d)               # calculate trajectory
    if d["WriteDataToFile"] == "yes":
        particle.WriteToFile(E_analytic, B_analytic, d)                                         # write data concerning the trajectory of the particle to a file if enabled in the xml-file
    return particle

pool = mp.Pool(processes=3)
results = [pool.apply_async(multiprocess, args=(particle, )) for particle in particles]
output = [p.get() for p in results]
print(output)



if len(particles) > 1:                  # if there are multiple particles, all of their trajectorys will also be plotted in a single figure
    PlotTrajectorys(output)"""









if len(particles) > 1:                  # if there are multiple particles, all of their trajectorys will also be plotted in a single figure
    PlotTrajectorys(particles)


if d["ElectricFieldPlot"] == "yes":                         # plot the electric field if enabled in the xml-file
    Plotfield(E_analytic, "Electric", d)

if d["MagneticFieldPlot"] == "yes":                         # plot the magnetic field if enabled in the xml-file
    Plotfield(B_analytic, "Magnetic", d)


input("Press Enter to end program and close all figures")








"""
Github:
# username: RobinBrabants
https://www.howtoforge.com/tutorial/install-git-and-github-on-ubuntu-14.04/
cd /home/robin/Documents/code
git add *
git commit -m "ee"
git push origin master
"""


# packages:
"""
python3.6 -m pip install scipy
python3.6 -m pip install sympy
python3.6 -m pip install matplotlib
apt-get install python3.6-tk
"""

# run as:
# python3.6 Main.py DataFile.xml
