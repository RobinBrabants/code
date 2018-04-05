from Lib.Functions import *
from sympy.vector import CoordSysCartesian, Vector
import sympy as sy


def ParticleMove(B, E, EWOS, d):

    print("ParticleMove")

    L = CoordSysCartesian('L')
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



    while time <= d["maxtime"]:

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

        Coordinates = [Coordinates[0] + vComponents[L.i]*d["t"], Coordinates[1] + vComponents[L.j]*d["t"], Coordinates[2] + vComponents[L.k]*d["t"]]
        Speed = (vComponents[L.i] + aComponents[L.i]*d["t"])*L.i + (vComponents[L.j] + aComponents[L.j]*d["t"])*L.j + (vComponents[L.k] + aComponents[L.k]*d["t"])*L.k
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