# Here all the electric and magnetic field producing elements will be represented by classes for which an analytical
# field can be calculated

import math
from abc import ABCMeta, abstractmethod
from sympy.vector import CoordSys3D, Vector
import math
import sympy as sy
from sympy.vector import CoordSys3D, Vector
from sympy import sqrt
from Lib.Functions import *
import numpy as np



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
    def GetClosestDistanceToPoint(self, point): # get closest distance from a point to this object
        pass

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
        CosBeta = ((x - x2) * (x2 - x1) + (y - y2) * (y2 - y1) + (z - z2) * (z2 - z1)) / (
        sy.sqrt((x - x2) ** 2 + (y - y2) ** 2 + (z - z2) ** 2) * sy.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2))
        n = ((y2 - y1) * (z - z1) - (z2 - z1) * (y - y1)) * L.i + ((z2 - z1) * (x - x1) - (x2 - x1) * (z - z1)) * L.j + ((x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)) * L.k

        eB = n.normalize()

        B = eB * (((mu * I) / (4 * pi * M)) * (CosAlpha - CosBeta))

        return B

    def GetClosestDistanceToPoint(self, point):
        return 0

    def IsPointInObject(self, point):
        return 0

    def FieldType(self):
        return 'magnetic'


#######################################################################################################################
#  electric field producing elements (potential on the surface us used):


class Sphere(Electrode):
    def __init__(self, name, Center, Radius, V):
        Electrode.__init__(self, name)
        self.center = Center
        self.radius = Radius
        self.potential = V

    def GetClosestDistanceToPoint(self, point):
        return abs(math.sqrt((self.center[0]-point[0])**2+(self.center[1]-point[1])**2+(self.center[2]-point[2])**2)-self.radius)

    def IsPointInObject(self, point):
        return (self.center[0]-point[0])**2+(self.center[1]-point[1])**2+(self.center[2]-point[2])**2 <= self.radius**2

    def FieldType(self):
        return 'electric'
