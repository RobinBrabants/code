from scipy.integrate import quad
import xml.etree.ElementTree as ET
from Lib.MagneticFieldComponents import *
from Lib.ElectricFieldComponents import *
from sympy import sqrt, diff, solve, nsolve
from scipy.optimize import fsolve, fmin
import numpy as np
import math
from random import random
from operator import itemgetter
import multiprocessing as mp
from Lib.Elements import *
from Lib.Objects_3D import *
from multiprocessing import pool
import time
import inspect
import sys


def ConvertAnglesToVector(Phi, Theta):
    # Function which converts given angles in radians into a normalised vector (see documentation for more info)
    pi = math.pi
    a = math.tan(Theta)/math.sqrt(1 + math.tan(Theta)**2)

    if 0 <= Phi <= pi / 2 or (3 / 2) * pi <= Phi < 2 * pi:
        ex = abs(math.cos(Phi) * a)
    else:
        ex = -abs(math.cos(Phi) * a)

    if 0 <= Phi <= pi:
        ey = abs(math.sqrt(a**2-ex**2))
    else:
        ey = -abs(math.sqrt(a ** 2 - ex ** 2))

    if math.tan(Theta) != 0:
        if 0 <= Theta <= pi/2:
            ez = abs(a/math.tan(Theta))
        else:
            ez = -abs(a / math.tan(Theta))
    else:
        ez = 1

    L = CoordSys3D('L')

    return ex*L.i + ey*L.j + ez*L.k


def getvector(length, Theta, Phi):
    # Function which converts given angles in radians into a vector with a certain length (see documentation for more info)
    eI = ConvertAnglesToVector(math.radians(Phi), math.radians(Theta))

    return eI*length


def ConvertVectorToAngles(v):
    # Function which converts a vector into angles in degrees (see documentation for more info)
    vcomponents = v.components
    UpdateDictionary(vcomponents)

    L = CoordSys3D('L')

    ex = vcomponents[L.i]
    ey = vcomponents[L.j]
    ez = vcomponents[L.k]

    a = math.sqrt(ey**2 + ex**2)
    Phi = math.degrees(math.acos(ex/a))
    Theta = math.degrees(math.atan(a/ez))

    return Phi, Theta


def UpdateDictionary(Dict):
    # Function which adds missing elements to the Dictionary object
    L = CoordSys3D('L')
    if L.i not in Dict:
        Dict.setdefault(L.i, 0)
    if L.j not in Dict:
        Dict.setdefault(L.j, 0)
    if L.k not in Dict:
        Dict.setdefault(L.k, 0)


def EvaluateAnalyticField(F, S):
    # Functions which evaluates an analytical field
    x, y, z, Phi2, t = sy.symbols('x y z Phi2 t')
    L = CoordSys3D('L')

    Fcomponents = F.components
    UpdateDictionary(Fcomponents)

    for basis in [L.i, L.j, L.k]:
        # For loop which checks whether there are any integrals which need to be numerically evaluated

        sum_integrals = 0

        if str(Fcomponents[basis]).find("Integral(") != -1:
            while str(str(Fcomponents[basis])).find("Integral(") != -1:
                str(F)
                # This somehow solved the problem that F would change throughout the process
                end1 = str(str(Fcomponents[basis])).find("Integral(")

                integrand = str(Fcomponents[basis])[str(Fcomponents[basis]).find("Integral(") + 9:str(Fcomponents[basis]).find(", ")]
                integrand = sy.sympify(integrand)
                integrand = integrand.subs([(x, S[0]), (y, S[1]), (z, S[2])])

                info1 = str(Fcomponents[basis])[str(Fcomponents[basis]).find(",") + 1:]
                dx = info1[info1.find("(") + 1:info1.find(",")]
                info2 = info1[info1.find(",") + 2:]
                a = info2[:info2.find(",")]
                b = info2[info2.find(",") + 2:info2.find("))")]

                begin2 = info2.find("))")


                if dx == "Phi2":
                    sum_integrals += quad(lambda Phi2: eval(str(integrand)), eval(a), eval(b))[0]
                if dx == "t":
                    sum_integrals += sy.integrate(eval(str(integrand)), (t,  eval(a), eval(b)))

                Fcomponents[basis] = str(Fcomponents[basis])[:end1] + info2[begin2 + 3:]

            if sum_integrals != 0:
                Fcomponents[basis] += str(sum_integrals)

            Fcomponents[basis] = eval(Fcomponents[basis])


    Fdic = Fcomponents[L.i]*L.i + Fcomponents[L.j]*L.j + Fcomponents[L.k]*L.k

    if Fdic != Vector.zero:
        Feval = Fdic.evalf(subs={x: S[0], y: S[1], z: S[2]})
    else:
        Feval = Vector.zero

    return Feval


def ReadXml():
    # User needs to put a named xml file on the same level as the Main.py executable (so in the program folder)
    # The name of the file needs to be passed when executing the program from the command window, e.g.: python Main.py Datafile.xml
    # output is a list of objects which represent the electrodes(_WOS), a particle object and a dictionary holding information about the setup

    # Check if the correct argument is given:
    if len(sys.argv) == 1:
        sys.exit("ERROR: User needs to give the xml-file holding the data as an argument")
    elif len(sys.argv) > 2:
        sys.exit("ERROR: Too many arguments given, please only give the data file as an argument")
    else:
        print("The data wil be read from: " + str(sys.argv[1]) + " ...")

    tree = ET.parse(str(sys.argv[1]))  # Load the given xml-file passed as an argument
    root = tree.getroot()

    def GetObjects(class_names, root):
        # function which will generate an object list from the enabled elements in the xml-file

        object_list = []

        for cls in class_names:
            if root.find('.//' + cls) is not None:  # checks if there are any elements in the xml-file which have the same name as the classes
                input_parameters_names = list(inspect.signature(eval(cls)).parameters)  # the input parameters needed as input to generate an object from a specific class
                input_parameters_names.remove('name')

                i = 1
                for element in root.iter(cls):  # if there are multiple elements from the same type, there will be iterated over these elements
                    if element.attrib["status"] == "enabled":  # will only make an object if the status is enabled in the xml file
                        input_parameters = [cls + str(i)]  # first input parameter is the object name
                        for parameter in input_parameters_names:
                            if element.find(parameter) is None:  # error message to warn for misspelling or forgetting to specify certain input parameters
                                sys.exit("ERROR: the input parameter " + str(parameter) + " for the element: " + str(cls) + " has not been specified or has been misspelled in the xml-file, please correct this and execute the program again")
                            input_parameters.append(eval(element.find(parameter).text))  # will find the other input parameters in the xml file and hold its values in a list

                        object_list.append(eval(cls)(*input_parameters))  # generates the correct class object with the extracted input parameters and appends it to the list
                        i += 1

        return object_list


    # ELECTRODES (electric and magnetic):
    class_names = [cls.__name__ for cls in vars()["Electrode"].__subclasses__()]  # get the class names derived from the electrode base class

    electrodes = GetObjects(class_names, root)  # list that holds all the specified electrodes as class objects

    # ELECTRODES_WOS (for the Walk on Spheres method, electric)
    class_names = [cls.__name__ + "_WOS" for cls in vars()["Object_3D"].__subclasses__()]  # get the class names derived from the Object_3D base class (and add _WOS)

    electrodes_WOS = GetObjects(class_names, root)  # list that holds all the specified electrodes_WOS as class objects

    # PARTICLE:
    class_name = "Particle"
    particle = 0

    """
    Particle = root.find("Particle")
    d["v"] = eval(Particle.find("Velocity").text)
    d["Theta1"] = eval(Particle.find("Theta1").text)
    d["Phi1"] = eval(Particle.find("Phi1").text)
    d["a"] = eval(Particle.find("Acceleration").text)
    d["Theta2"] = eval(Particle.find("Theta2").text)
    d["Phi2"] = eval(Particle.find("Phi2").text)

    Position = Particle.find("Position")
    d["Position"] = [eval(Position.find("x").text), eval(Position.find("y").text), eval(Position.find("z").text)]

    if Particle.find("Type").text == "custom":
        d["Mass"] = eval(Particle.find("Mass").text)
        d["Charge"] = eval(Particle.find("Charge").text)
    else:
        if Particle.find("Type").text in ParticleDictionary:
            d["Mass"] = ParticleDictionary[Particle.find("Type").text][0]
            d["Charge"] = ParticleDictionary[Particle.find("Type").text][1]
        else:
            print("Particle not found in dictionary")
    """

    # SETUP (works with dictionary instead of class objects)
    d = {}          # empty dictionary to which names will be appended which will later be used in different classes/functions

    Setup = root.find("Setup")

    Trajectory = Setup.find("Trajectory")
    TrajectoryBoundaries = Trajectory.find("TrajectoryBoundaries")
    d["xmin"] = eval(TrajectoryBoundaries.find("xmin").text)
    d["xmax"] = eval(TrajectoryBoundaries.find("xmax").text)
    d["ymin"] = eval(TrajectoryBoundaries.find("ymin").text)
    d["ymax"] = eval(TrajectoryBoundaries.find("ymax").text)
    d["zmin"] = eval(TrajectoryBoundaries.find("zmin").text)
    d["zmax"] = eval(TrajectoryBoundaries.find("zmax").text)
    d["timesteps"] = eval(Trajectory.find("TimeSteps").text)
    d["timelimit"] = eval(Trajectory.find("TimeLimit").text)
    d["interval"] = eval(Trajectory.find("Interval").text)

    Output = Setup.find("Output")
    d["TrajectoryPlot"] = Output.find("TrajectoryPlot").attrib["execute"]
    d["MagneticFieldPlot"] = Output.find("MagneticFieldPlot").attrib["execute"]
    d["NormalizeMagneticFieldPlot"] = Output.find("MagneticFieldPlot").attrib["normalize"]
    d["ElectricFieldPlot"] = Output.find("ElectricFieldPlot").attrib["execute"]
    d["NormalizeElectricFieldPlot"] = Output.find("ElectricFieldPlot").attrib["normalize"]
    d["WriteDataToFile"] = Output.find("WriteDataToFile").attrib["execute"]

    MagneticFieldPlot = Output.find("MagneticFieldPlot")
    MagneticFieldBoundaries = MagneticFieldPlot.find("MagneticFieldBoundaries")
    d["xmin1"] = eval(MagneticFieldBoundaries.find("xmin").text)
    d["xmax1"] = eval(MagneticFieldBoundaries.find("xmax").text)
    d["ymin1"] = eval(MagneticFieldBoundaries.find("ymin").text)
    d["ymax1"] = eval(MagneticFieldBoundaries.find("ymax").text)
    d["zmin1"] = eval(MagneticFieldBoundaries.find("zmin").text)
    d["zmax1"] = eval(MagneticFieldBoundaries.find("zmax").text)

    ElectricFieldPlot = Output.find("ElectricFieldPlot")
    ElectricFieldBoundaries = ElectricFieldPlot.find("ElectricFieldBoundaries")
    d["xmin2"] = eval(ElectricFieldBoundaries.find("xmin").text)
    d["xmax2"] = eval(ElectricFieldBoundaries.find("xmax").text)
    d["ymin2"] = eval(ElectricFieldBoundaries.find("ymin").text)
    d["ymax2"] = eval(ElectricFieldBoundaries.find("ymax").text)
    d["zmin2"] = eval(ElectricFieldBoundaries.find("zmin").text)
    d["zmax2"] = eval(ElectricFieldBoundaries.find("zmax").text)

    d["FileName"] = Output.find("WriteDataToFile").text







    # ask if given setup is correct before proceding







    return electrodes, electrodes_WOS, particle, d




def ResultingField(electrodes):
    # function which iterates over the objects initialised by the ReadXml function in order to get an analytic expression for the magnetic and electric vectorfield

    B = Vector.zero
    E = Vector.zero

    for electrode in electrodes:
        if electrode.FieldType() =="magnetic":
            B += electrode.GetField()
        elif electrode.FieldType() =="electric":
            E += electrode.GetField()

    return B, E
