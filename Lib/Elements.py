# Here all the electric and magnetic field producing elements will be represented by classes for which an analytical
# field can be calculated


import math
from abc import abstractmethod
from sympy.vector import CoordSys3D, Vector
import sympy as sy
from Lib.Objects_3D import *


class Electrode(object):
    # a class to represent all electrodes from which to derive the actual shaped objects later
    def __init__(self, name):
        self.name = name                        # give up a name for the object

    def PrintName(self):                        # output name and some info on the object
        print('name: ', self.name);             # everything needs to be implemented

    @abstractmethod
    def GetField(self):                         # get the magnetic or electric field vector as a function of the
        pass                                    # 3D coordinate points x, y and z

    @abstractmethod
    def IsPointInObject(self, point):           # see whether a point is inside an object
        pass

    @abstractmethod
    def FieldType(self):                        # see what kind of field electrodes produce (magnetic or electric)
        pass



# now instantiate specific classes and give them their specific functionality by overwriting the generic function


#######################################################################################################################
#  magnetic field producing elements (watch out: conventional current directions are used):


class StraightConductor(Electrode):
    # a straight line conductor
    def __init__(self, name, CoordinatesPoint1, CoordinatesPoint2, Current, Radius):
        Electrode.__init__(self, name)
        self.CoordinatesPoint1 = CoordinatesPoint1      # starting point straight conductor
        self.CoordinatesPoint2 = CoordinatesPoint2      # ending point straight conductor
        self.Current = Current                          # current
        self.Radius = Radius                            # radius of the straightconductor (represented as a cylinder)

        # Conventional current direction goes from point 1 to point 2, unless otherwise specified (see example xml-file)

    def GetField(self):
        # see the documentation on how the field is calculated (Magnetisch veld eindige rechte geleider)
        P = self.CoordinatesPoint1
        Q = self.CoordinatesPoint2
        I = self.Current

        x1 = P[0]
        y1 = P[1]
        z1 = P[2]
        x2 = Q[0]
        y2 = Q[1]
        z2 = Q[2]

        pi = math.pi
        mu = 4 * pi * 10 ** (-7)

        L = CoordSys3D('L')

        x, y, z = sy.symbols('x y z')

        x3 = x1 + ((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1) + (z2 - z1) * (z - z1)) / ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) * (x2 - x1)
        y3 = y1 + ((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1) + (z2 - z1) * (z - z1)) / ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) * (y2 - y1)
        z3 = z1 + ((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1) + (z2 - z1) * (z - z1)) / ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) * (z2 - z1)

        M = sy.sqrt((x - x3) ** 2 + (y - y3) ** 2 + (z - z3) ** 2)

        CosAlpha = ((x - x1) * (x2 - x1) + (y - y1) * (y2 - y1) + (z - z1) * (z2 - z1)) / (sy.sqrt((x - x1) ** 2 + (y - y1) ** 2 + (z - z1) ** 2) * sy.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2))
        CosBeta = ((x - x2) * (x2 - x1) + (y - y2) * (y2 - y1) + (z - z2) * (z2 - z1)) / (sy.sqrt((x - x2) ** 2 + (y - y2) ** 2 + (z - z2) ** 2) * sy.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2))
        n = ((y2 - y1) * (z - z1) - (z2 - z1) * (y - y1)) * L.i + ((z2 - z1) * (x - x1) - (x2 - x1) * (z - z1)) * L.j + ((x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)) * L.k

        eB = n.normalize()

        B = eB * (((mu * I) / (4 * pi * M)) * (CosAlpha - CosBeta))

        return B

    def IsPointInObject(self, point, interval):
        return 0

    def FieldType(self):
        return 'magnetic'


class StraightConductorCollection(Electrode):
    # an assembly of straight line conductors
    def __init__(self, name, Current, Radius, Points):
        Electrode.__init__(self, name)
        self.Current = Current                          # current
        self.Radius = Radius                            # radius of the straightconductors (represented as a cylinders)
        self.Points = Points                            # starting and end points of the seperate straight line conductors
                                                        # these points will all be connected by straight conductors to form a collection of straight conductors

        # Conventional current direction goes from point 1 to point 2, unless otherwise specified (see example xml-file)

    def GetField(self):
        B = Vector.zero

        # initialize separate straight conductors and add their vector fields to get the total magnetic field:
        for i in range(0, len(self.Points) - 1):
            straightconductor = StraightConductor("name", self.Points[i], self.Points[i+1], self.Current, self.Radius)
            B += straightconductor.GetField()
        return B

    def IsPointInObject(self, point, interval):
        return 0

    def FieldType(self):
        return 'magnetic'


class RectangularCoil(Electrode):
    # a rectangular coil existing of straight line conductors
    def __init__(self, name, Current, Radius, Phi, Theta, Psi, Length, Width, Heigth, StartingPoint, Windings, ClockWise):
        Electrode.__init__(self, name)
        self.Current = Current                          # current
        self.Radius = Radius                            # radius of the straightconductors (represented as a cylinders)
        self.Phi = Phi                                  # these are the three Euler Angles to define the orientation of the coil (in degrees) (see documentation)
        self.Theta = Theta
        self.Psi = Psi
        self.Length = Length                            # dimension coil in e_1 direction
        self.Width = Width                              # dimension coil in e_2 direction
        self.Heigth = Heigth                            # dimension coil in e_3 direction
        self.StartingPoint = StartingPoint              # origin of the axial system of the coil
        self.Windings = Windings                        # number of windings (is counted in quarters of windings, e.g. 21.75 windings would be a valid argument)
        self.ClockWise = ClockWise                      # whether the windings of the coil should start clockwise (argument = True --> clockwise, argument = False --> anti-clockwise)

        # Conventional current direction goes from point 1 to point 2, unless otherwise specified (see example xml-file)

    def GetField(self):
        # see the documentation on how the field is calculated (Magnetisch veld rechthoekige spoel)
        from Lib.Functions import UpdateDictionary
        from math import sin, cos

        L = CoordSys3D('L')

        Phi = math.radians(self.Phi)
        Theta = math.radians(self.Theta)
        Psi = math.radians(self.Psi)

        L = CoordSys3D('L')

        e_1 = (cos(Phi) * cos(Psi) - sin(Phi) * sin(Psi) * cos(Theta)) * L.i + (sin(Phi) * cos(Psi) + cos(Phi) * sin(Psi) * cos(Theta)) * L.j + sin(Theta) * sin(Psi) * L.k
        e_2 = (-cos(Phi) * sin(Psi) - sin(Phi) * cos(Psi) * cos(Theta)) * L.i + (-sin(Phi) * sin(Psi) + cos(Phi) * cos(Psi) * cos(Theta)) * L.j + sin(Theta) * cos(Psi) * L.k
        e_3 = (sin(Phi) * sin(Theta)) * L.i + (-cos(Phi) * sin(Theta)) * L.j + cos(Theta) * L.k


        if self.ClockWise == "True":
            a_1 = self.Width * e_1
            a_2 = self.Length * e_2
            a_3_1 = ((self.Width / (2 * self.Width + 2 * self.Length)) * (self.Heigth / self.Windings)) * e_3
            a_3_2 = ((self.Length / (2 * self.Width + 2 * self.Length)) * (self.Heigth / self.Windings)) * e_3
        else:
            a_1 = self.Length * e_1
            a_2 = self.Width * e_2
            a_3_1 = ((self.Length / (2 * self.Width + 2 * self.Length)) * (self.Heigth / self.Windings)) * e_3
            a_3_2 = ((self.Width / (2 * self.Width + 2 * self.Length)) * (self.Heigth / self.Windings)) * e_3

        a_1 = a_1.components
        UpdateDictionary(a_1)
        a_2 = a_2.components
        UpdateDictionary(a_2)
        a_3_1 = a_3_1.components
        UpdateDictionary(a_3_1)
        a_3_2 = a_3_2.components
        UpdateDictionary(a_3_2)

        new_point = self.StartingPoint
        windings = 0
        list = []
        list.append(self.StartingPoint[:])

        while windings < self.Windings:
            new_point[0] += a_1[L.i] + a_3_1[L.i]
            new_point[1] += a_1[L.j] + a_3_1[L.j]
            new_point[2] += a_1[L.k] + a_3_1[L.k]
            list.append(new_point[:])
            windings += 0.25

            if windings >= self.Windings:
                break

            new_point[0] += a_2[L.i] + a_3_2[L.i]
            new_point[1] += a_2[L.j] + a_3_2[L.j]
            new_point[2] += a_2[L.k] + a_3_2[L.k]
            list.append(new_point[:])
            windings += 0.25

            if windings >= self.Windings:
                break

            new_point[0] += -a_1[L.i] + a_3_1[L.i]
            new_point[1] += -a_1[L.j] + a_3_1[L.j]
            new_point[2] += -a_1[L.k] + a_3_1[L.k]
            list.append(new_point[:])
            windings += 0.25

            if windings >= self.Windings:
                break

            new_point[0] += -a_2[L.i] + a_3_2[L.i]
            new_point[1] += -a_2[L.j] + a_3_2[L.j]
            new_point[2] += -a_2[L.k] + a_3_2[L.k]
            list.append(new_point[:])
            windings += 0.25

        collection = StraightConductorCollection("name", self.Current, self.Radius, list)


        B = collection.GetField()

        return B

    def IsPointInObject(self, point, interval):
        return 0

    def FieldType(self):
        return 'magnetic'



#######################################################################################################################
#  electric field producing elements (potential on the surface is used):


class Sphere_Field(Electrode):
    # a sphere on which a potential is applied
    def __init__(self, name, CoordinatesCenter, Radius, Potential):
        Electrode.__init__(self, name)
        self.CoordinatesCenter = CoordinatesCenter      # coordinates of the center of a sphere
        self.Radius = Radius                            # radius of the sphere
        self.Potential = Potential                      # potential on the surface of the sphere

    def GetField(self):
        # law of Coulomb is used to calculate the electrical field

        pi = math.pi
        e0 = 8.854187817 * 10 ** (-12)      # permitivity of vacuum

        M = self.CoordinatesCenter
        x1 = M[0]
        y1 = M[1]
        z1 = M[2]

        L = CoordSys3D('L')

        x, y, z = sy.symbols('x y z')

        distance = sy.sqrt((x - x1) ** 2 + (y - y1) ** 2 + (z - z1) ** 2)      # distance squared between center of sphere
                                                                            # and point in space
        r = (x - x1) * L.i + (y - y1) * L.j + (z - z1) * L.k
        er = r.normalize()                                      # gives correct direction of electrical field

        Q = self.Potential*4*pi*e0*distance                     # charge on sphere calcuted from potential on sphere
                                                                # (V = Q/(4.pi.e0.distance))

        E = Q / (4 * pi * e0 * distance**2) * er           # law of Coulomb

        return E

    def IsPointInObject(self, point, interval):

        sphere = Sphere(self.name, self.CoordinatesCenter, self.Radius, self.Potential)

        return sphere.IsPointInObject(point, interval)

    def FieldType(self):
        return 'electric'
