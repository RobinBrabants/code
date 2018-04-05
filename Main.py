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


#d = ReadXml()

electrodes = []

tree = ET.parse("DataFile.xml")
root = tree.getroot()

d = {}


class_names = [cls.__name__ for cls in vars()['Electrode'].__subclasses__()]



for cls in class_names:
    if root.find('.//' + cls):
        print (cls)

        for element in root.iter(cls):
            if element.attrib["status"] == "enabled":
                print("ok")
                print(StraightConductor.__dict__.iteritems())


"""
for straightconductors in root.iter("StraightConductor"):
    if straightconductors.attrib["status"] == "enabled":

        CoordinatesPoint = straightconductors.find("CoordinatesPoint1")
        CoordinatesPoint1 = [eval(CoordinatesPoint.find("x").text), eval(CoordinatesPoint.find("y").text),
                            eval(CoordinatesPoint.find("z").text)]

        CoordinatesPoint = straightconductors.find("CoordinatesPoint2")
        CoordinatesPoint2 = [eval(CoordinatesPoint.find("x").text), eval(CoordinatesPoint.find("y").text),
                            eval(CoordinatesPoint.find("z").text)]

        Current = eval(straightconductors.find("Current").text)
        Radius = eval(straightconductors.find("Radius").text)

        print(CoordinatesPoint1, CoordinatesPoint2, Current, Radius)

        Element = StraightConductor("StraightConductor_" + str(i), CoordinatesPoint1, CoordinatesPoint2, Current, Radius)
        electrodes.append(Element)
        i += 1

for electrode in electrodes:
    print(electrode.name)
    print(electrode.GetClosestDistanceToPoint([0,0,11]))
    print(electrode.Current)
"""
"""
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
