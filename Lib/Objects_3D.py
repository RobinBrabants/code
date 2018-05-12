# Here all the 3 dimensional objects which are being used will be stored as classes with the main goal of calculating
# the minimal distances to these objects and to determine whether a point is inside an object

# classes from Elements will only use the IsPointInObject function to determine in the ParticleMove function from the
# Particle class whether a particle had a collision with an object.

# The DistanceToObject function will be used for the Walk On Spheres method. This is explained in the
# functions_WOS file.


import math
from math import sqrt
from abc import abstractmethod
from sympy.vector import CoordSys3D


class Object_3D(object):
    # a class to represent all the objects from which to derive the actual shaped objects later
    def __init__(self, name):
        self.name = name                        # give up a name for the object

    @abstractmethod
    def DistanceToObject(self, point): # get closest distance from a point to this object (negative distance
                                                # means the point is inside the object), used in function IsPointInObject
        pass

    def IsPointInObject(self, point, interval): # see whether a point is inside an object and determines whether a given
                                                # point in space is near enough to the surface of an object to be
                                                # considered as colliding with the object

        # 0 means particle is not near enough to object to be considered colliding with it and is not inside object
        # 1 means particle is considered to be colliding with the object (near enough)
        # 2 means particle is inside object

        distance = self.DistanceToObject(point)

        if distance > interval:
            return 0
        elif 0 <= distance <= interval:
            return 1
        else:
            return 2



# now instantiate specific classes and give them their specific functionality by overwriting the generic function:

#####################################################################################################################


class Sphere(Object_3D):
    # represents a 3D sphere surface
    def __init__(self, name, CoordinatesCenter, Radius, Potential):
        Object_3D.__init__(self, name)
        self.CoordinatesCenter = CoordinatesCenter      # coordinates of the center of a sphere
        self.Radius = Radius                            # radius of the sphere
        self.Potential = Potential                      # potential on the surface of the sphere


    def DistanceToObject(self, point):
        # self explanatory calculation of the minimal distance to a sphere
        # determines distance between a given point in space and the center of the sphere and then substracts the radius

        M = self.CoordinatesCenter
        x_1 = M[0]
        y_1 = M[1]
        z_1 = M[2]

        x = point[0]
        y = point[1]
        z = point[2]

        distance = abs(sqrt((x_1 - x) ** 2 + (y_1 - y) ** 2 + (z_1 - z) ** 2)) - self.Radius

        return distance


class CircularDisk(Object_3D):
    # represents a 3D infinitely small circular disk
    def __init__(self, name, CoordinatesCenter, Radius, Theta, Phi, Potential):
        Object_3D.__init__(self, name)
        self.CoordinatesCenter = CoordinatesCenter      # coordinates of the center of the disk
        self.Radius = Radius                            # radius of the disk
        self.Theta = Theta                              # measured from the x-axis in anti-clockwise direction until the projection of the normal circle vector on the xy-plane (0-360 degree)
        self.Phi = Phi                                  # measured from the projection of the velocity vector on the xy-plane in anti-clockwise direction until the normal circle vector (0-360 degree)
        self.Potential = Potential                      # potential on the surface of the disk


    def DistanceToObject(self, point):
        # see the documentation on how the distance is calculated (Afstand tot een cirkelschijf)
        from Lib.Functions import ConvertAnglesToVector, UpdateDictionary

        Phi = math.radians(self.Phi)
        Theta = math.radians(self.Theta)
        N = ConvertAnglesToVector(Theta, Phi)           #  normal circle vector

        M = self.CoordinatesCenter
        x_1 = M[0]
        y_1 = M[1]
        z_1 = M[2]

        x = point[0]
        y = point[1]
        z = point[2]

        L = CoordSys3D('L')

        delta = (x - x_1) * L.i + (y - y_1) * L.j + (z - z_1) * L.k

        P = (x) * L.i + (y) * L.j + (z) * L.k
        C =  (x_1) * L.i + (y_1) * L.j + (z_1) * L.k
        Q = delta - (N.dot(delta))*N + C
        D = Q-C

        if D.magnitude() <= self.Radius:
            distance = (P - Q).magnitude()
        else:
            distance = sqrt((N.dot(delta)) ** 2 + (abs((N.cross(delta)).magnitude()) - self.Radius) ** 2)

        return distance


class Cylinder(Object_3D):
    # represents a 3D cylinder closed at the bottom and top
    def __init__(self, name, CoordinatesPoint1, CoordinatesPoint2, Radius, Potential):
        Object_3D.__init__(self, name)
        self.CoordinatesPoint1 = CoordinatesPoint1      # starting point cylinder
        self.CoordinatesPoint2 = CoordinatesPoint2      # ending point cylinder
        self.Radius = Radius                            # radius of the cylinder
        self.Potential = Potential                      # potential on the surface of the cylinder


    def DistanceToObject(self, point):
        from Lib.Functions import UpdateDictionary

        P = self.CoordinatesPoint1
        Q = self.CoordinatesPoint2

        x_1 = P[0]
        y_1 = P[1]
        z_1 = P[2]

        x_2 = Q[0]
        y_2 = Q[1]
        z_2 = Q[2]

        x = point[0]
        y = point[1]
        z = point[2]

        L = CoordSys3D('L')

        v = (x_2 - x_1) * L.i + (y_2 - y_1) * L.j + (z_2 - z_1) * L.k
        e_v = v.normalize()
        e_v_c = e_v.components
        UpdateDictionary(e_v_c)
        a = e_v_c[L.i]
        b = e_v_c[L.j]
        c = e_v_c[L.k]

        if c != 0:
            if z_1 > z_2:
                z_1 = z_2
                z_2 = z_1
            if (1 / c) * (a * x_1 + b * y_1 + c * z_1 - a * x - b * y) <= z <= (1 / c) * (a * x_2 + b * y_2 + c * z_2 - a * x - b * y):
                t = (a * x + b * y + c * z - a * x_1 - b * y_1 - z_1) / (
                c * (z_2 - z_1) + a * (x_2 - x_1) + b * (y_2 - y_1))
                x_c = x_1 + t * (x_2 - x_1)
                y_c = y_1 + t * (y_2 - y_1)
                z_c = z_1 + t * (z_2 - z_1)
                distance = abs(sqrt((x - x_c) ** 2 + (y - y_c) ** 2 + (z - z_c) ** 2)) - self.Radius
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
                    distance = abs(sqrt((x - x_c) ** 2 + (y - y_c) ** 2 + (z - z_c) ** 2)) - self.Radius
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
                    distance = abs(sqrt((x - x_c) ** 2 + (y - y_c) ** 2 + (z - z_c) ** 2)) - self.Radius
                elif x > (1 / a) * (a * x_2 + b * y_2 + c * z_2 - b * y - c * z):
                    delta = (x - x_2) * L.i + (y - y_2) * L.j + (z - z_2) * L.k
                    distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - self.Radius) ** 2)
                else:
                    delta = (x - x_1) * L.i + (y - y_1) * L.j + (z - z_1) * L.k
                    distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - self.Radius) ** 2)

        return distance
