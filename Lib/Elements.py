# Here all the magnetic field producing elements will be represented by classes for which an analytical
# field can be calculated


import math
from abc import abstractmethod
from sympy.vector import CoordSys3D, Vector
import sympy as sy
from Lib.Objects_3D import *


class FieldSource(object):
    # a class to represent all electrodes from which to derive the actual shaped objects later
    def __init__(self, name):
        self.name = name                        # give up a name for the object

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


class StraightConductor(FieldSource):
    # a straight line conductor
    def __init__(self, name, CoordinatesPoint1, CoordinatesPoint2, Current, Radius):
        FieldSource.__init__(self, name)
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


class StraightConductorCollection(FieldSource):
    # (a.k.a. piecewise linear)s
    # an assembly of straight line conductors
    def __init__(self, name, Current, Radius, Points):
        FieldSource.__init__(self, name)
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


class RectangularCoil(FieldSource):
    # a rectangular coil existing of straight line conductors
    def __init__(self, name, Current, Radius, Phi, Theta, Psi, Length, Width, Heigth, StartingPoint, Windings, ClockWise):
        FieldSource.__init__(self, name)
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


class BentConductor(FieldSource):
    # a bent conductor
    def __init__(self, name, CoordinatesCenter, CurvatureRadius, Theta, Phi, Interval, Current, Radius):
        FieldSource.__init__(self, name)
        self.CoordinatesCenter = CoordinatesCenter      # coordinates of the center of the bent conductor
        self.CurvatureRadius = CurvatureRadius          # radius of curvature
        self.Theta = Theta                              # measured from the x-axis in anti-clockwise direction until the projection of the normal circle vector on the xy-plane (0-360 degree)
        self.Phi = Phi                                  # measured from the projection of the velocity vector on the xy-plane in anti-clockwise direction until the normal circle vector (0-360 degree)
        self.Interval = Interval                        # angle starts counting from a perpendicular vector on the normal circle vector parallel to the xy-plane
        self.Current = Current                          # current
        self.Radius = Radius                            # radius, thickness of the bent conductor

        # The current in the conductor flows as if the normal circle vector were its rotation vector, to reverse make the Current negative

    def GetField(self):
        # see the documentation on how the field is calculated (Magnetisch veld cirkelvormige geleider)

        pi = math.pi
        mu = 4 * pi * 10 ** (-7)

        M = self.CoordinatesCenter
        R = self.CurvatureRadius

        x1 = M[0]
        y1 = M[1]
        z1 = M[2]

        Theta = pi/2 - math.radians(self.Phi)
        Phi = math.radians(self.Theta)

        L = CoordSys3D('L')

        x, y, z, Phi2 = sy.symbols('x y z Phi2')

        dl = R * (sy.sin(Phi2) * sy.sin(Phi) + sy.cos(Phi2) * sy.cos(Theta) * sy.cos(Phi)) * L.i + R * (-sy.sin(Phi2) * sy.cos(Phi) + sy.cos(Phi2) * sy.cos(Theta) * sy.sin(Phi)) * L.j + R * (-sy.cos(Phi2) * sy.sin(Theta)) * L.k

        r = (x - x1 + R * sy.cos(Phi2) * sy.sin(Phi) - R * sy.sin(Phi2) * sy.cos(Theta) * sy.cos(Phi)) * L.i + (y - y1 - R * sy.cos(Phi2) * sy.cos(Phi) - R * sy.sin(Phi2) * sy.cos(Theta) * sy.sin(Phi)) * L.j + (z - z1 + R * sy.sin(Phi2) * sy.sin(Theta)) * L.k

        b = ((mu * self.Current) / (4 * pi)) * ((dl.cross(r)) / (r.magnitude() ** 3))

        B = sy.integrate(-b, (Phi2, math.radians(self.Interval[0]), math.radians(self.Interval[1])))

        return B

    def IsPointInObject(self, point, interval):
        return 0

    def FieldType(self):
        return 'magnetic'


class CircularConductor(FieldSource):
    # a circular conductor
    def __init__(self, name, CoordinatesCenter, CurvatureRadius, Theta, Phi, Current, Radius):
        FieldSource.__init__(self, name)
        self.CoordinatesCenter = CoordinatesCenter      # coordinates of the center of the bent conductor
        self.CurvatureRadius = CurvatureRadius          # radius of curvature
        self.Theta = Theta                              # measured from the x-axis in anti-clockwise direction until the projection of the normal circle vector on the xy-plane (0-360 degree)
        self.Phi = Phi                                  # measured from the projection of the velocity vector on the xy-plane in anti-clockwise direction until the normal circle vector (0-360 degree)
        self.Current = Current                          # current
        self.Radius = Radius                            # radius, thickness of the bent conductor

        # The current in the conductor flows as if the normal circle vector were its rotation vector, to reverse make the Current negative

    def GetField(self):
        # see the documentation on how the field is calculated (Magnetisch veld cirkelvormige geleider)

        circle = BentConductor("name", self.CoordinatesCenter, self.CurvatureRadius, self.Theta, self.Phi, [0, 360], self.Current, self.Radius)

        return circle.GetField()

    def IsPointInObject(self, point, interval):
        return 0

    def FieldType(self):
        return 'magnetic'


class CircularCoil(FieldSource):
    # a circular conductor
    def __init__(self, name, CoordinatesCenter, CurvatureRadius, Theta, Phi, Heigth, Windings, BeginAngle, Current, Radius):
        FieldSource.__init__(self, name)
        self.CoordinatesCenter = CoordinatesCenter      # coordinates of the center of the bottom of the coil
        self.CurvatureRadius = CurvatureRadius          # radius of curvature of the coil itself, not its windings
        self.Theta = Theta                              # measured from the x-axis in anti-clockwise direction until the projection of the coil orientation vector on the xy-plane (0-180 degree)
        self.Phi = Phi                                  # measured from the projection of the velocity vector on the xy-plane in anti-clockwise direction until the coil orientation vector (0-180 degree)
        self.Heigth = Heigth                            # heigth of the coil
        self.Windings = Windings                        # number of windings (is counted in halves of windings, e.g. 21.5 windings would be a valid argument)
        self.BeginAngle = BeginAngle                    # angle starts counting from a perpendicular vector on the normal circle vector parallel to the xy-plane, this is for the first winding
        self.Current = Current                          # current
        self.Radius = Radius                            # radius, thickness of the bent conductors

        # Conventional current direction goes from point 1 to point 2, unless otherwise specified (see example xml-file)

    def GetField(self):
        # see the documentation on how the field is calculated (Magnetisch veld cirkelvormige geleider)

        from Lib.Functions import UpdateDictionary, ConvertAnglesToVector, ConvertVectorToAngles
        from sympy import solve_poly_system

        Theta = math.radians(self.Theta)
        Phi = math.radians(self.Phi)
        M = self.CoordinatesCenter

        pi = math.pi
        mu = 4 * pi * 10 ** (-7)

        h = 0.25 * (self.Heigth / self.Windings)

        L = CoordSys3D('L')

        n = ConvertAnglesToVector(Theta, Phi)

        M = M[0] * L.i + M[1] * L.j + M[2] * L.k
        M = M + h * n

        Mcomponents = M.components
        UpdateDictionary(Mcomponents)

        M = [Mcomponents[L.i], Mcomponents[L.j], Mcomponents[L.k]]

        ncomponents = n.components
        UpdateDictionary(ncomponents)

        nx = ncomponents[L.i]
        ny = ncomponents[L.j]
        nz = ncomponents[L.k]

        kx, ky = sy.symbols('kx ky')

        if nx != 0 and ny != 0:
            k = solve_poly_system([kx * nx + ky * ny, kx ** 2 + ky ** 2 - 1], kx, ky)
            kx = k[0][0]
            ky = k[0][1]
        else:
            kx = 0
            ky = 1

        tan_alpha = h / self.CurvatureRadius

        dnx, dny, dnz = sy.symbols('dnx, dny, dnz')

        dn = solve_poly_system([nx * dnx + ny * dny + nz * dnz, dnx ** 2 + dny ** 2 + dnz ** 2 - tan_alpha ** 2,(kx / tan_alpha) * dnx + (ky / tan_alpha) * dny - math.cos(self.BeginAngle)], dnx, dny, dnz)

        if len(dn) >= 2:
            dnx = dn[1][0]
            dny = dn[1][1]
            dnz = dn[1][2]
        else:
            dnx = dn[0][0]
            dny = dn[0][1]
            dnz = dn[0][2]

        dn = dnx * L.i + dny * L.j + dnz * L.k

        r = n + dn

        r = r.normalize()

        Phi, Theta = ConvertVectorToAngles(r)

        alpha = math.atan(tan_alpha)

        R = self.CurvatureRadius / math.cos(alpha)

        B = Vector.zero
        windings = 0

        begin = self.BeginAngle

        while windings < self.Windings:
            bentconductor = BentConductor("name", M, R, Theta, Phi, [begin, begin + 180], self.Current, self.Radius)
            B += bentconductor.GetField()

            M = M[0] * L.i + M[1] * L.j + M[2] * L.k
            M = M + 2 * h * n
            Mcomponents = M.components
            UpdateDictionary(Mcomponents)
            M = [Mcomponents[L.i], Mcomponents[L.j], Mcomponents[L.k]]

            dn = -dn
            r = n + dn
            r = r.normalize()
            Theta, Phi = ConvertVectorToAngles(r)

            begin = begin + 180

            windings += 0.5

        return B

    def IsPointInObject(self, point, interval):
        return 0

    def FieldType(self):
        return 'magnetic'



#######################################################################################################################
#  electric field producing elements (potential on the surface is used):
# This is just a test function for the Walk On Spheres method


class Sphere_Field(FieldSource):
    # a sphere on which a potential is applied
    def __init__(self, name, CoordinatesCenter, Radius, Potential):
        FieldSource.__init__(self, name)
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
