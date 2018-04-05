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
from multiprocessing import pool
import time


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

    L = CoordSysCartesian('L')

    return ex*L.i + ey*L.j + ez*L.k


def getvector(length, Theta, Phi):
    # Function which converts given angles in radians into a vector with a certain length (see documentation for more info)
    eI = ConvertAnglesToVector(math.radians(Phi), math.radians(Theta))

    return eI*length


def ConvertVectorToAngles(v):
    # Function which converts a vector into angles in degrees (see documentation for more info)
    vcomponents = v.components
    UpdateDictionary(vcomponents)

    L = CoordSysCartesian('L')

    ex = vcomponents[L.i]
    ey = vcomponents[L.j]
    ez = vcomponents[L.k]

    a = math.sqrt(ey**2 + ex**2)
    Phi = math.degrees(math.acos(ex/a))
    Theta = math.degrees(math.atan(a/ez))

    return Phi, Theta


def UpdateDictionary(Dict):
    # Function which adds missing elements to the Dictionary object
    L = CoordSysCartesian('L')
    if L.i not in Dict:
        Dict.setdefault(L.i, 0)
    if L.j not in Dict:
        Dict.setdefault(L.j, 0)
    if L.k not in Dict:
        Dict.setdefault(L.k, 0)


def EvaluateAnalyticField(F, S):
    # Functions which evaluates an analytical field
    x, y, z, Phi2, t = sy.symbols('x y z Phi2 t')
    L = CoordSysCartesian('L')

    Fcomponents = F.components
    UpdateDictionary(Fcomponents)

    for basis in [L.i, L.j, L.k]:
        # For loop which checks whether there are any integrals that need to be numerically evaluated

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
    #funtion which reads the XML file
    tree = ET.parse("DataFile.xml")
    root = tree.getroot()

    d = {}
    ParticleDictionary = {"electron": [9.10938356*(10**(-31)), -1.6021766208*(10**(-19))], "proton": [1.672621898*(10**(-27)), 1.6021766208*10**(-19)]}
    i = 1

    """
    # function to get the data (list of electrode classes) for the main function
    electrodes = []

    i = 1
    while "Msphere2_" + str(i) in d:
        electrodes.append(
            Sphere("sphere_" + str(i), d["Msphere2_" + str(i)], d["Rsphere2_" + str(i)], d["Vsphere2_" + str(i)]))
        i += 1

    i = 1
    while "Pcylinder_" + str(i) in d:
        electrodes.append(Cylinder(d["Pcylinder_" + str(i)], d["Qcylinder_" + str(i)], d["Rcylinder_" + str(i)],
                                   d["Vcylinder_" + str(i)]))
        i += 1

    i = 1
    while "Mcircdisc_" + str(i) in d:
        electrodes.append(
            CircularDisk(d["Mcircdisc_" + str(i)], d["Rcircdisc_" + str(i)], d["Phicircdisc_" + str(i)],
                         d["Thetacircdisc_" + str(i)], d["Vcircdisc_" + str(i)]))
        i += 1

    i = 1
    while "Mcircdisc4h_" + str(i) in d:
        electrodes.append(
            CircularDisk4Holes(d["Mcircdisc4h_" + str(i)], d["Mcircdisc4h1_" + str(i)], d["Mcircdisc4h2_" + str(i)],
                               d["Mcircdisc4h3_" + str(i)], d["Mcircdisc4h4_" + str(i)], d["Rcircdisc4h_" + str(i)],
                               d["Rcircdisc4h1_" + str(i)], d["Rcircdisc4h2_" + str(i)],
                               d["Rcircdisc4h3_" + str(i)], d["Rcircdisc4h4_" + str(i)], d["Ncircdisc4h_" + str(i)],
                               d["Vcircdisc4h_" + str(i)]))
        i += 1

    """


    """
    for StraightConductor in root.iter("StraightConductor"):
        if StraightConductor.attrib["status"] == "enabled":
            CoordinatesPoint1 = StraightConductor.find("CoordinatesPoint1")
            d["P_" + str(i)] = [eval(CoordinatesPoint1.find("x").text), eval(CoordinatesPoint1.find("y").text), eval(CoordinatesPoint1.find("z").text)]
            CoordinatesPoint2 = StraightConductor.find("CoordinatesPoint2")
            d["Q_" + str(i)] = [eval(CoordinatesPoint2.find("x").text), eval(CoordinatesPoint2.find("y").text), eval(CoordinatesPoint2.find("z").text)]
            d["I_" + str(i)] = eval(StraightConductor.find("Current").text)
            i += 1


    for StraightConductorCollection in root.iter("StraightConductorCollection"):
        if StraightConductorCollection.attrib["status"] == "enabled":
            amount = StraightConductorCollection.attrib["amount"]
            for j in range(1, int(amount) + 1):
                CoordinatesPoint = StraightConductorCollection.find("CoordinatesPoint" +str(j))
                d[str(i) + "S_" + str(j)] = [eval(CoordinatesPoint.find("x").text), eval(CoordinatesPoint.find("y").text),
                                eval(CoordinatesPoint.find("z").text)]
            d["Icollection_" + str(i)] = eval(StraightConductorCollection.find("Current").text)
            i += 1


    for RectangularCoil in root.iter("RectangularCoil"):
        if RectangularCoil.attrib["status"] == "enabled":
            StartingPoint = RectangularCoil.find("StartingPoint")
            d["w_" + str(i)] = eval(RectangularCoil.find("Width").text)
            d["l_" + str(i)] = eval(RectangularCoil.find("Length").text)
            d["h_" + str(i)] = eval(RectangularCoil.find("Heigth").text)
            d["N_" + str(i)] = eval(RectangularCoil.find("Windings").text)
            d["SP_" + str(i)] = [eval(StartingPoint.find("x").text), eval(StartingPoint.find("y").text),
                                 eval(StartingPoint.find("z").text)]
            d["Phireccoil_" + str(i)] = eval(RectangularCoil.find("Phi").text)
            d["Thetareccoil_" + str(i)] = eval(RectangularCoil.find("Theta").text)
            d["Psireccoil_" + str(i)] = eval(RectangularCoil.find("Psi").text)
            d["Ireccoil_" + str(i)] = eval(RectangularCoil.find("Current").text)
            d["begin_" + str(i)] = eval(RectangularCoil.find("begin").text)
            i += 1

    for CircularConductor in root.iter("CircularConductor"):
        if CircularConductor.attrib["status"] == "enabled":
            CoordinatesCentre = CircularConductor.find("CoordinatesCentre")
            d["M_" + str(i)] = [eval(CoordinatesCentre.find("x").text), eval(CoordinatesCentre.find("y").text), eval(CoordinatesCentre.find("z").text)]
            d["R_" + str(i)] = eval(CircularConductor.find("Radius").text)
            Orientation = CircularConductor.find("Orientation")
            d["Phi_" + str(i)] = math.radians(eval(Orientation.find("Phi").text))
            d["Theta_" + str(i)] = math.radians(eval(Orientation.find("Theta").text))
            d["Icircle_" + str(i)] = eval(CircularConductor.find("Current").text)
            i += 1



    for BentConductor in root.iter("BentConductor"):
        if BentConductor.attrib["status"] == "enabled":
            CoordinatesCentre = BentConductor.find("CoordinatesCentre")
            Interval = BentConductor.find("Interval")
            d["Mbent_" + str(i)] = [eval(CoordinatesCentre.find("x").text), eval(CoordinatesCentre.find("y").text), eval(CoordinatesCentre.find("z").text)]
            d["Rbent_" + str(i)] = eval(BentConductor.find("Radius").text)
            Orientation = BentConductor.find("Orientation")
            d["Phibent_" + str(i)] = math.radians(eval(Orientation.find("Phi").text))
            d["Thetabent_" + str(i)] = math.radians(eval(Orientation.find("Theta").text))
            d["Ibent_" + str(i)] = eval(BentConductor.find("Current").text)
            d["Interval_" + str(i)] = [eval(Interval.find("from").text), eval(Interval.find("to").text)]
            i += 1



    for CircularCoil in root.iter("CircularCoil"):
        if CircularCoil.attrib["status"] == "enabled":
            CoordinatesCentreCoil = CircularCoil.find("CoordinatesCentre")
            d["hcoil_" + str(i)] = eval(CircularCoil.find("Heigth").text)
            d["Ncoil_" + str(i)] = eval(CircularCoil.find("Windings").text)
            d["Mcoil_" + str(i)] = [eval(CoordinatesCentreCoil.find("x").text), eval(CoordinatesCentreCoil.find("y").text),
                                 eval(CoordinatesCentreCoil.find("z").text)]
            d["Rcoil_" + str(i)] = eval(CircularCoil.find("Radius").text)
            d["Phicoil_" + str(i)] = eval(CircularCoil.find("Phi").text)
            d["Thetacoil_" + str(i)] = eval(CircularCoil.find("Theta").text)
            d["Icoil_" + str(i)] = eval(CircularCoil.find("Current").text)
            d["begincoil_" + str(i)] = eval(CircularCoil.find("begin").text)
            i += 1


    for Sphere in root.iter("Sphere"):
        if Sphere.attrib["status"] == "enabled":
            CoordinatesCentre = Sphere.find("CoordinatesCentre")
            d["Msphere_" + str(i)] = [eval(CoordinatesCentre.find("x").text), eval(CoordinatesCentre.find("y").text), eval(CoordinatesCentre.find("z").text)]
            d["Rsphere_" + str(i)] = eval(Sphere.find("Radius").text)
            d["Qsphere_" + str(i)] = eval(Sphere.find("Charge").text)
            i += 1


    for Sphere2 in root.iter("Sphere2"):
        if Sphere2.attrib["status"] == "enabled":
            CoordinatesCentre = Sphere2.find("CoordinatesCentre")
            d["Msphere2_" + str(i)] = [eval(CoordinatesCentre.find("x").text), eval(CoordinatesCentre.find("y").text), eval(CoordinatesCentre.find("z").text)]
            d["Rsphere2_" + str(i)] = eval(Sphere2.find("Radius").text)
            d["Vsphere2_" + str(i)] = eval(Sphere2.find("Potential").text)

            i += 1

    for Cylinder in root.iter("Cylinder"):
        if Cylinder.attrib["status"] == "enabled":
            CoordinatesPoint1 = Cylinder.find("CoordinatesPoint1")
            d["Pcylinder_" + str(i)] = [eval(CoordinatesPoint1.find("x").text), eval(CoordinatesPoint1.find("y").text), eval(CoordinatesPoint1.find("z").text)]
            CoordinatesPoint2 = Cylinder.find("CoordinatesPoint2")
            d["Qcylinder_" + str(i)] = [eval(CoordinatesPoint2.find("x").text), eval(CoordinatesPoint2.find("y").text), eval(CoordinatesPoint2.find("z").text)]
            d["Rcylinder_" + str(i)] = eval(Cylinder.find("Radius").text)
            d["Vcylinder_" + str(i)] = eval(Cylinder.find("Potential").text)
            i += 1


    for CircularDisc in root.iter("CircularDisc"):
        if CircularDisc.attrib["status"] == "enabled":
            CoordinatesCentre = CircularDisc.find("CoordinatesCentre")
            d["Mcircdisc_" + str(i)] = [eval(CoordinatesCentre.find("x").text), eval(CoordinatesCentre.find("y").text), eval(CoordinatesCentre.find("z").text)]
            d["Rcircdisc_" + str(i)] = eval(CircularDisc.find("Radius").text)
            Orientation = CircularDisc.find("Orientation")
            d["Phicircdisc_" + str(i)] = math.radians(eval(Orientation.find("Phi").text))
            d["Thetacircdisc_" + str(i)] = math.radians(eval(Orientation.find("Theta").text))
            d["Vcircdisc_" + str(i)] = eval(CircularDisc.find("Potential").text)

            i += 1


    for CircularDisc4Holes in root.iter("CircularDisc4Holes"):
        if CircularDisc4Holes.attrib["status"] == "enabled":
            CoordinatesCentre = CircularDisc4Holes.find("CoordinatesCentre")
            d["Mcircdisc4h_" + str(i)] = [eval(CoordinatesCentre.find("x").text), eval(CoordinatesCentre.find("y").text), eval(CoordinatesCentre.find("z").text)]
            d["Rcircdisc4h_" + str(i)] = eval(CircularDisc4Holes.find("Radius").text)
            Orientation = CircularDisc4Holes.find("Orientation")
            d["Ncircdisc4h_" + str(i)] = [eval(Orientation.find("x").text), eval(Orientation.find("y").text), eval(Orientation.find("z").text)]
            d["Vcircdisc4h_" + str(i)] = eval(CircularDisc4Holes.find("Potential").text)
            Holes = CircularDisc4Holes.find("Holes")
            d["Rcircdisc4h1_" + str(i)] = eval(Holes.find("Radius1").text)
            d["Rcircdisc4h2_" + str(i)] = eval(Holes.find("Radius2").text)
            d["Rcircdisc4h3_" + str(i)] = eval(Holes.find("Radius3").text)
            d["Rcircdisc4h4_" + str(i)] = eval(Holes.find("Radius4").text)
            CoordinatesCentre1 = Holes.find("CoordinatesCentre1")
            d["Mcircdisc4h1_" + str(i)] = [eval(CoordinatesCentre1.find("x").text),eval(CoordinatesCentre1.find("y").text),eval(CoordinatesCentre1.find("z").text)]
            CoordinatesCentre2 = Holes.find("CoordinatesCentre2")
            d["Mcircdisc4h2_" + str(i)] = [eval(CoordinatesCentre2.find("x").text),eval(CoordinatesCentre2.find("y").text),eval(CoordinatesCentre2.find("z").text)]
            CoordinatesCentre3 = Holes.find("CoordinatesCentre3")
            d["Mcircdisc4h3_" + str(i)] = [eval(CoordinatesCentre3.find("x").text),eval(CoordinatesCentre3.find("y").text),eval(CoordinatesCentre3.find("z").text)]
            CoordinatesCentre4 = Holes.find("CoordinatesCentre4")
            d["Mcircdisc4h4_" + str(i)] = [eval(CoordinatesCentre4.find("x").text),eval(CoordinatesCentre4.find("y").text),eval(CoordinatesCentre4.find("z").text)]

            i += 1


    for Surface in root.iter("Surface"):
        if Surface.attrib["status"] == "enabled":
            x, y = sy.symbols('x y')
            d["SurfEq_" + str(i)] = eval(Surface.find("Equation").text)
            d["Vsurface_" + str(i)] = eval(Surface.find("Potential").text)
            i += 1


    for Line in root.iter("Line"):
        if Line.attrib["status"] == "enabled":
            CoordinatesPoint1 = Line.find("CoordinatesPoint1")
            d["Pline_" + str(i)] = [eval(CoordinatesPoint1.find("x").text), eval(CoordinatesPoint1.find("y").text),
                                eval(CoordinatesPoint1.find("z").text)]
            CoordinatesPoint2 = Line.find("CoordinatesPoint2")
            d["Rline_" + str(i)] = [eval(CoordinatesPoint2.find("x").text), eval(CoordinatesPoint2.find("y").text),
                                eval(CoordinatesPoint2.find("z").text)]
            d["Qline_" + str(i)] = eval(Line.find("Charge").text)
            i += 1

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


    Setup = root.find("Setup")

    Trajectory = Setup.find("Trajectory")
    TrajectoryBoundaries = Trajectory.find("TrajectoryBoundaries")
    d["xmin"] = eval(TrajectoryBoundaries.find("xmin").text)
    d["xmax"] = eval(TrajectoryBoundaries.find("xmax").text)
    d["ymin"] = eval(TrajectoryBoundaries.find("ymin").text)
    d["ymax"] = eval(TrajectoryBoundaries.find("ymax").text)
    d["zmin"] = eval(TrajectoryBoundaries.find("zmin").text)
    d["zmax"] = eval(TrajectoryBoundaries.find("zmax").text)
    d["t"] = eval(Trajectory.find("TimeSteps").text)
    d["maxtime"] = eval(Trajectory.find("TimeLimit").text)

    Output = Setup.find("Output")
    d["TrajectoryPlot"] = Output.find("TrajectoryPlot").attrib["execute"]
    d["MagneticFieldPlot"] = Output.find("MagneticFieldPlot").attrib["execute"]
    d["NormalizeMagneticFieldPlot"] = Output.find("MagneticFieldPlot").attrib["normalize"]
    d["ElectricFieldPlot"] = Output.find("ElectricFieldPlot").attrib["execute"]
    d["NormalizeElectricFieldPlot"] = Output.find("ElectricFieldPlot").attrib["normalize"]
    d["WriteDataToFile"] = Output.find("WriteDataToFile").attrib["execute"]

    MagneticFieldPlot = Output.find("MagneticFieldPlot")
    MagneticFieldBoundaries =  MagneticFieldPlot.find("MagneticFieldBoundaries")
    d["xmin1"] = eval(MagneticFieldBoundaries.find("xmin").text)
    d["xmax1"] = eval(MagneticFieldBoundaries.find("xmax").text)
    d["ymin1"] = eval(MagneticFieldBoundaries.find("ymin").text)
    d["ymax1"] = eval(MagneticFieldBoundaries.find("ymax").text)
    d["zmin1"] = eval(MagneticFieldBoundaries.find("zmin").text)
    d["zmax1"] = eval(MagneticFieldBoundaries.find("zmax").text)

    ElectricFieldPlot = Output.find("ElectricFieldPlot")
    ElectricFieldBoundaries =  ElectricFieldPlot.find("ElectricFieldBoundaries")
    d["xmin2"] = eval(ElectricFieldBoundaries.find("xmin").text)
    d["xmax2"] = eval(ElectricFieldBoundaries.find("xmax").text)
    d["ymin2"] = eval(ElectricFieldBoundaries.find("ymin").text)
    d["ymax2"] = eval(ElectricFieldBoundaries.find("ymax").text)
    d["zmin2"] = eval(ElectricFieldBoundaries.find("zmin").text)
    d["zmax2"] = eval(ElectricFieldBoundaries.find("zmax").text)

    d["FileName"] = Output.find("WriteDataToFile").text


    return d


def SetupElements(d):
    i = 1
    Blist = []


    while "P_" + str(i) in d:
        Blist.append(StraightConductor(d["P_" + str(i)], d["Q_" + str(i)], d["I_" + str(i)]))
        i += 1


    i = 1
    j = 1

    while str(j) + "S_" + str(i) in d:
        list = ()
        while str(j) + "S_" + str(i) in d:
            list += (d[str(j) + "S_" + str(i)],)
            i += 1
        Blist.append(StraightConductorCollection(d["Icollection_" + str(j)], *list))
        i = 1
        j += 1


    i = 1
    while "SP_" + str(i) in d:
        Blist.append(Rectangularcoil(d["w_" + str(i)], d["l_" + str(i)], d["h_" + str(i)], d["N_" + str(i)], d["SP_" + str(i)], d["Phireccoil_" + str(i)],  d["Thetareccoil_" + str(i)], d["Psireccoil_" + str(i)], d["Ireccoil_" + str(i)], d["begin_" + str(i)]))
        i += 1

    i = 1

    while "M_" + str(i) in d:
        Blist.append(CircularConductor(d["M_" + str(i)], d["R_" + str(i)], d["Phi_" + str(i)], d["Theta_" + str(i)], d["Icircle_" + str(i)]))
        i += 1


    i = 1

    while "Mbent_" + str(i) in d:
        Blist.append(BentConductor(d["Mbent_" + str(i)], d["Rbent_" + str(i)], d["Phibent_" + str(i)], d["Thetabent_" + str(i)], d["Interval_" + str(i)], d["Ibent_" + str(i)]))
        i += 1


    i = 1

    while "Mcoil_" + str(i) in d:
        Blist.append(CircularCoil(d["Mcoil_" + str(i)], d["Rcoil_" + str(i)], d["Phicoil_" + str(i)], d["Thetacoil_" + str(i)], d["begincoil_" + str(i)], d["hcoil_" + str(i)], d["Ncoil_" + str(i)],  d["Icoil_" + str(i)]))
        i += 1


    return Blist


def SetupElements2(d):
    i = 1
    Elist = []


    while "Msphere_" + str(i) in d:
        Elist.append(Sphere(d["Msphere_" + str(i)], d["Rsphere_" + str(i)], d["Qsphere_" + str(i)]))
        i += 1


    i = 1
    while "Pline_" + str(i) in d:
        Elist.append(Line(d["Pline_" + str(i)], d["Rline_" + str(i)], d["Qline_" + str(i)]))
        i += 1


    return Elist


def ResultingField(list):
    F = Vector.zero

    for value in list:
       F += value

    return F

def ResultingString(list):
    F = ""
    if list != None:

        for value in list:
           F += value
    return F
