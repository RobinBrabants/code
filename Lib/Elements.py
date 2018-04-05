# Here all the electric and magnetic field producing elements will be represented by classes for which an analytical
# field can be calculated

import math
from abc import ABCMeta, abstractmethod
from sympy.vector import CoordSysCartesian, Vector
import math
import sympy as sy
from sympy.vector import CoordSysCartesian, Vector
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
#  magnetic field producing elements:

class StraightConductor(Electrode):
    def __init__(self, name, CoordinatesPoint1, CoordinatesPoint2, Current, Radius):
        Electrode.__init__(self, name)
        self.CoordinatesPoint1 = CoordinatesPoint1
        self.CoordinatesPoint2 = CoordinatesPoint2
        self.Current = Current
        self.Radius = Radius

    def GetField(self):
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

        L = CoordSysCartesian('L')

        x, y, z = sy.symbols('x y z')

        x3 = x1 + ((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1) + (z2 - z1) * (z - z1)) / (
        (x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) * (x2 - x1)
        y3 = y1 + ((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1) + (z2 - z1) * (z - z1)) / (
        (x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) * (y2 - y1)
        z3 = z1 + ((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1) + (z2 - z1) * (z - z1)) / (
        (x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) * (z2 - z1)

        M = sy.sqrt((x - x3) ** 2 + (y - y3) ** 2 + (z - z3) ** 2)

        CosAlpha = ((x - x1) * (x2 - x1) + (y - y1) * (y2 - y1) + (z - z1) * (z2 - z1)) / (
        sy.sqrt((x - x1) ** 2 + (y - y1) ** 2 + (z - z1) ** 2) * sy.sqrt(
            (x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2))
        CosBeta = ((x - x2) * (x2 - x1) + (y - y2) * (y2 - y1) + (z - z2) * (z2 - z1)) / (
        sy.sqrt((x - x2) ** 2 + (y - y2) ** 2 + (z - z2) ** 2) * sy.sqrt(
            (x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2))
        n = ((y2 - y1) * (z - z1) - (z2 - z1) * (y - y1)) * L.i + (
                                                                  (z2 - z1) * (x - x1) - (x2 - x1) * (z - z1)) * L.j + (
                                                                                                                       (
                                                                                                                       x2 - x1) * (
                                                                                                                       y - y1) - (
                                                                                                                       y2 - y1) * (
                                                                                                                       x - x1)) * L.k

        eB = n.normalize()

        B = eB * (((mu * I) / (4 * pi * M)) * (CosAlpha - CosBeta))

        return B

    def GetClosestDistanceToPoint(self, point):
        from Lib.Functions import UpdateDictionary
        from math import sqrt

        x = point[0]
        y = point[1]
        z = point[2]
        x_1 = self.CoordinatesPoint1[0]
        y_1 = self.CoordinatesPoint1[1]
        z_1 = self.CoordinatesPoint1[2]
        x_2 = self.CoordinatesPoint2[0]
        y_2 = self.CoordinatesPoint2[1]
        z_2 = self.CoordinatesPoint2[2]

        L = CoordSysCartesian('L')

        v = (x_2 - x_1) * L.i + (y_2 - y_1) * L.j + (z_2 - z_1) * L.k
        e_v = v.normalize()
        e_v_c = e_v.components
        UpdateDictionary(e_v_c)
        a = e_v_c[L.i]
        b = e_v_c[L.j]
        c = e_v_c[L.k]

        if c != 0:
            if (1 / c) * (a * x_1 + b * y_1 + c * z_1 - a * x - b * y) <= z <= (1 / c) * (
                                a * x_2 + b * y_2 + c * z_2 - a * x - b * y):
                t = (a * x + b * y + c * z - a * x_1 - b * y_1 - z_1) / (
                c * (z_2 - z_1) + a * (x_2 - x_1) + b * (y_2 - y_1))
                x_c = x_1 + t * (x_2 - x_1)
                y_c = y_1 + t * (y_2 - y_1)
                z_c = z_1 + t * (z_2 - z_1)
                distance = abs(sqrt((x - x_c) ** 2 + (y - y_c) ** 2 + (z - z_c) ** 2) - self.Radius)
            elif z > (1 / c) * (a * x_2 + b * y_2 + c * z_2 - a * x - b * y):
                delta = (x - x_2) * L.i + (y - y_2) * L.j + (z - z_2) * L.k
                distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - self.Radius) ** 2)
            else:
                delta = (x - x_1) * L.i + (y - y_1) * L.j + (z - z_1) * L.k
                distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - self.Radius) ** 2)
        else:
            if b != 0:
                if (1 / b) * (a * x_1 + b * y_1 + c * z_1 - a * x - c * z) <= y <= (1 / b) * (
                                    a * x_2 + b * y_2 + c * z_2 - a * x - c * z):
                    t = (a * x + b * y + c * z - a * x_1 - b * y_1 - z_1) / (
                    c * (z_2 - z_1) + a * (x_2 - x_1) + b * (y_2 - y_1))
                    x_c = x_1 + t * (x_2 - x_1)
                    y_c = y_1 + t * (y_2 - y_1)
                    z_c = z_1 + t * (z_2 - z_1)
                    distance = abs(sqrt((x - x_c) ** 2 + (y - y_c) ** 2 + (z - z_c) ** 2) - self.Radius)
                elif y > (1 / b) * (a * x_2 + b * y_2 + c * z_2 - a * x - c * z):
                    delta = (x - x_2) * L.i + (y - y_2) * L.j + (z - z_2) * L.k
                    distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - self.Radius) ** 2)
                else:
                    delta = (x - x_1) * L.i + (y - y_1) * L.j + (z - z_1) * L.k
                    distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - self.Radius) ** 2)
            else:
                if (1 / a) * (a * x_1 + b * y_1 + c * z_1 - b * y - c * z) <= x <= (1 / a) * (
                                    a * x_2 + b * y_2 + c * z_2 - b * y - c * z):
                    t = (a * x + b * y + c * z - a * x_1 - b * y_1 - z_1) / (
                    c * (z_2 - z_1) + a * (x_2 - x_1) + b * (y_2 - y_1))
                    x_c = x_1 + t * (x_2 - x_1)
                    y_c = y_1 + t * (y_2 - y_1)
                    z_c = z_1 + t * (z_2 - z_1)
                    distance = abs(sqrt((x - x_c) ** 2 + (y - y_c) ** 2 + (z - z_c) ** 2) - self.Radius)
                elif x > (1 / a) * (a * x_2 + b * y_2 + c * z_2 - b * y - c * z):
                    delta = (x - x_2) * L.i + (y - y_2) * L.j + (z - z_2) * L.k
                    distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - self.Radius) ** 2)
                else:
                    delta = (x - x_1) * L.i + (y - y_1) * L.j + (z - z_1) * L.k
                    distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - self.Radius) ** 2)

        return distance

    def IsPointInObject(self, point):
        return (self.center[0]-point[0])**2+(self.center[1]-point[1])**2+(self.center[2]-point[2])**2 <= self.radius**2

    def FieldType(self):
        return 'magnetic'


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
