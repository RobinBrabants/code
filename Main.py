# main function from which everything in the Lib will be called and the program will be executed


from Lib.Particle import *
from Lib.RepresentData import *
from Lib.Functions import *
import sys




"""from Lib.Functions_WOS import *

Point = [0, 0, 10]

electrodes, electrodes_WOS, particles, d = ReadXml()        # extract all the data from the xml file
print(electrodes)
print(electrodes_WOS)

B_analytic, E_analytic = ResultingField(electrodes)
E, B = GetFields(Point, Vector.zero, B_analytic, E_analytic, electrodes_WOS)

print("\r\nanalytical: " + str(E))
print("WOS: " + str(ElectricalField_WOS(electrodes_WOS, Point)))"""



electrodes, electrodes_WOS, particles, d = ReadXml()        # extract all the data from the xml file

for particle in particles:
    # check whether the starting position of the particle lies out of the predetermined box
    if particle.Position[0] < d["xmin"] or particle.Position[0] > d["xmax"] or particle.Position[1] < d["ymin"] or particle.Position[1] > d["ymax"] or particle.Position[2] < d["zmin"] or particle.Position[2] > d["zmax"]:
        sys.exit("ERROR: The Starting position of %s (%s) lies out of the predetermined box, check your xml-file\r\n" % (particle.name, particle.Type))

    # check if the starting position of the particle is inside an object
    for electrode in electrodes + electrodes_WOS:
        if electrode.IsPointInObject(particle.Position, d["interval"]) == 2:
            sys.exit("ERROR: The Starting position of %s (%s) lies inside %s, check your xml-file\n" % (particle.name, particle.Type, electrode.name))


B_analytic, E_analytic = ResultingField(electrodes)         # already calculates the analytical fields for those objects for which it is possible

"""for particle in particles:
    trajectory = particle.ParticleMove(B_analytic, E_analytic, electrodes, electrodes_WOS, d)               # calculate trajectory
    if d["WriteDataToFile"] == "yes":
        particle.WriteToFile(E_analytic, B_analytic, trajectory, d)                                         # write data concerning the trajectory of the particle to a file if enabled in the xml-file
    if d["TrajectoryPlot"] == "yes":
        particle.PlotTrajectory(trajectory)"""                                           # plot the trajectory if enabled in the xml-file

if d["ElectricFieldPlot"] == "yes":                         # plot the electric field if enabled in the xml-file
    Plotfield(E_analytic, "Electric", d)

if d["MagneticFieldPlot"] == "yes":                         # plot the magnetic field if enabled in the xml-file
    Plotfield(B_analytic, "Magnetic", d)

input("Press Enter to end program and close all figures")











"""
https://www.howtoforge.com/tutorial/install-git-and-github-on-ubuntu-14.04/
cd /home/robin/Documents/code
git add *
git commit -m "ee"
git push origin master
"""
# username: RobinBrabants

# http://chris35wills.github.io/conda_python_version/

# https://pythonhosted.org/PyInstaller/usage.html


# packages:
"""
python3.6 -m pip install scipy
python3.6 -m pip install sympy
python3.6 -m pip install matplotlib
apt-get install python3.6-tk
"""

# run as:
# python3.6 Main.py DataFile.xml





















"""from Lib.Objects_3D import *

CoordinatesCenter = [0, 0, 0]
Radius = 5
Potential =4
sphere = Sphere("test", CoordinatesCenter, Radius, Potential)
point = [6.1, 0, 0]
print(sphere.GetClosestDistanceToPoint(point))
print(sphere.IsPointInObject(point,1))"""