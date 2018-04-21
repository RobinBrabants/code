from Lib.Functions import *
from sympy.vector import CoordSys3D, Vector
import sympy as sy

class Particle():        # represents a 3D sphere surface
    def __init__(self, name, CoordinatesCenter, Radius, Potential):
        electrodes, electrodes_WOS, particle, d = ReadXml()  # extract all the data from the xml file

        B_analytic, E_analytic = ResultingField(
            electrodes)  # already calculates the analytical fields for those objects for which it is possible

        self.name = name
        self.CoordinatesCenter = CoordinatesCenter      # coordinates of the center of a sphere
        self.Radius = Radius                            # radius of the sphere
        self.Potential = Potential                      # potential on the surface of the sphere


    def GetClosestDistanceToPoint(self, point):
        # self explanatory calculation of the minimal distance to a sphere
        # determines distance between a given point in space and the center of the sphere and then substracts the radius
        from math import sqrt

        M = self.CoordinatesCenter
        x_1 = M[0]
        y_1 = M[1]
        z_1 = M[2]

        x = point[0]
        y = point[1]
        z = point[2]

        distance = abs(sqrt((x_1 - x) ** 2 + (y_1 - y) ** 2 + (z_1 - z) ** 2) - self.Radius)

        return distance


def ParticleMove(B, E, EWOS, d):

    print("ParticleMove")

    L = CoordSys3D('L')
    x, y, z = sy.symbols('x y z')

    Coordinates = d["Position"]
    Speed = getvector(d["v"], d["Phi1"], d["Theta1"],)
    Accleration = getvector(d["a"], d["Phi2"], d["Theta2"])
    time = 0

    Trajectory = {"x": [], "y": [], "z": [], "t": [], "|v|": [], "|a|": []}
    Trajectory["x"].append(Coordinates[0])
    Trajectory["y"].append(Coordinates[1])
    Trajectory["z"].append(Coordinates[2])
    Trajectory["t"].append(time)
    Trajectory["|v|"].append(Speed.magnitude())
    Trajectory["|a|"].append(Accleration.magnitude())



    while time <= d["timelimit"]:

        B = EvaluateAnalyticField(B, Coordinates)
        E_analytic = EvaluateAnalyticField(E, Coordinates)
        #E_WOS =
        E = E_analytic



        #relativistic correction
        c = 299792458
        gamma = 1/(1-Speed.magnitude()**2/c**2)

        Blocal = gamma*(B - (Speed.cross(E))/c**2) - (gamma - 1)*(B.dot(Speed))*Speed.normalize()
        Elocal = gamma*(E + Speed.cross(B)) - (gamma - 1)*(E.dot(Speed))*Speed.normalize()

        vComponents = Speed.components

        # oftewel veld in 2 punten evalueren --> gemiddelde oftwel veld in punt tussen begin en eindpunt evalueren

        Accleration = (d["Charge"] / d["Mass"]) * ((Speed.cross(Blocal)) + Elocal) + Accleration
        aComponents = Accleration.components

        UpdateDictionary(vComponents)
        UpdateDictionary(aComponents)

        Coordinates = [Coordinates[0] + vComponents[L.i]*d["timesteps"], Coordinates[1] + vComponents[L.j]*d["timesteps"], Coordinates[2] + vComponents[L.k]*d["timesteps"]]
        Speed = (vComponents[L.i] + aComponents[L.i]*d["timesteps"])*L.i + (vComponents[L.j] + aComponents[L.j]*d["timesteps"])*L.j + (vComponents[L.k] + aComponents[L.k]*d["timesteps"])*L.k
        time += d["t"]

        if Coordinates[0] < d["xmin"] or Coordinates[0] > d["xmax"] or Coordinates[1] < d["ymin"] or Coordinates[1] > d["ymax"] or Coordinates[2] < d["zmin"] or Coordinates[2] > d["zmax"]:
            break

        Trajectory["x"].append(Coordinates[0])
        Trajectory["y"].append(Coordinates[1])
        Trajectory["z"].append(Coordinates[2])
        Trajectory["t"].append(time)
        Trajectory["|v|"].append(Speed.magnitude())
        Trajectory["|a|"].append(Accleration.magnitude())

    return Trajectory