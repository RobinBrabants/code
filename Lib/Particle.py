from Lib.Functions import *
from sympy.vector import CoordSys3D, Vector
import sympy as sy


ParticleDictionary = {"electron": [9.10938356 * (10 ** (-31)), -1.6021766208 * (10 ** (-19))], "proton": [1.672621898 * (10 ** (-27)), 1.6021766208 * 10 ** (-19)]}
# this is a dictionary which contains information about certain types of particles (can be expanded)
# in the form of: "Name": [mass, charge]

class Particle():
    def __init__(self, name, Type, Mass, Charge, Velocity, Theta1, Phi1, Acceleration, Theta2, Phi2, Position):
        self.name = name
        self.Type = Type        # the particle type (like electron, proton, muon...)

        if self.Type in ParticleDictionary:                  # If the particle is found in the ParticleDictionary the Mass and Charge of the particle will be correctly assigned
            self.Mass = ParticleDictionary[self.Type][0]
            self.Charge = ParticleDictionary[self.Type][1]
        else:
            self.Mass = Mass            # in kilograms
            self.Charge = Charge        # in coulombs

        self.Velocity = Velocity                # starting velocity of the particle in m/s
        self.Theta1 = Theta1                    # measured from the x-axis in anti-clockwise direction untill the projection of the velocity vector on the xy-plane (0-360 degree)
        self.Phi1 = Phi1                        # measured from the projection of the velocity vector on the xy-plane in anti-clockwise direction untill the velocity vector (0-360 degree)
        self.Acceleration = Acceleration        # acceleration of the particle which is the result of an external field and will all the time be added to the acceleration caused by the electric and magnetic fields (in m/s²)
        self.Theta2 = Theta2                    # measured from the x-axis in anti-clockwise direction untill the projection of the acceleration vector on the xy-plane (0-360 degree)
        self.Phi2 = Phi2                        # measured from the projection of the velocity vector on the xy-plane in anti-clockwise direction untill the acceleration vector (0-360 degree)
        self.Position = Position                # starting position particle

    def ParticleMove(self, B_analytic, E_analytic, electrodes, electrodes_WOS, d):
        print("The trajectory of the particle is now being calculated")

        L = CoordSys3D('L')
        x, y, z = sy.symbols('x y z')

        Coordinates = d["Position"]
        Speed = getvector(d["v"], d["Theta1"],  d["Phi1"])
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
            # E_WOS =
            E = E_analytic

            # relativistic correction
            c = 299792458
            gamma = 1 / (1 - Speed.magnitude() ** 2 / c ** 2)

            Blocal = gamma * (B - (Speed.cross(E)) / c ** 2) - (gamma - 1) * (B.dot(Speed)) * Speed.normalize()
            Elocal = gamma * (E + Speed.cross(B)) - (gamma - 1) * (E.dot(Speed)) * Speed.normalize()

            vComponents = Speed.components

            # oftewel veld in 2 punten evalueren --> gemiddelde oftwel veld in punt tussen begin en eindpunt evalueren

            Accleration = (d["Charge"] / d["Mass"]) * ((Speed.cross(Blocal)) + Elocal) + Accleration
            aComponents = Accleration.components

            UpdateDictionary(vComponents)
            UpdateDictionary(aComponents)

            Coordinates = [Coordinates[0] + vComponents[L.i] * d["timesteps"],
                           Coordinates[1] + vComponents[L.j] * d["timesteps"],
                           Coordinates[2] + vComponents[L.k] * d["timesteps"]]
            Speed = (vComponents[L.i] + aComponents[L.i] * d["timesteps"]) * L.i + (
                        vComponents[L.j] + aComponents[L.j] * d["timesteps"]) * L.j + (
                                vComponents[L.k] + aComponents[L.k] * d["timesteps"]) * L.k
            time += d["t"]

            if Coordinates[0] < d["xmin"] or Coordinates[0] > d["xmax"] or Coordinates[1] < d["ymin"] or Coordinates[
                1] > d["ymax"] or Coordinates[2] < d["zmin"] or Coordinates[2] > d["zmax"]:
                break

            Trajectory["x"].append(Coordinates[0])
            Trajectory["y"].append(Coordinates[1])
            Trajectory["z"].append(Coordinates[2])
            Trajectory["t"].append(time)
            Trajectory["|v|"].append(Speed.magnitude())
            Trajectory["|a|"].append(Accleration.magnitude())


        print("The calculations for the trajectory of the particle have been successfully finished")

        return Trajectory

    def PlotTrajectory(data):

        print("PlotTrajectory")

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        x = data["x"]
        y = data["y"]
        z = data["z"]
        t = data["t"]

        xs = np.array(x)
        ys = np.array(y)
        zs = np.array(z)
        ts = np.array(t)

        ax.plot(x, y, z)

        len_t = len(t) - 1
        ind_pos = []

        for i in range(0, 11):
            ind_pos.append(0 + i * (int(len_t / 10)))

        xx = (xs[ind_pos])
        yy = (ys[ind_pos])
        zz = (zs[ind_pos])
        tt = (ts[ind_pos])

        for t, x, y, z in zip(tt, xx, yy, zz):
            label = '%s' % t
            ax.text(x, y, z, label)
            ax.scatter(x, y, z, c="red")

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        fig.suptitle('trajectory of particle', fontsize=14, fontweight='bold')
        ax.set_title('time in seconds', style='italic', fontsize=8, bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 7})

        plt.show(block="False")

    def WriteFile(B, data, d):

        print("WriteFile")

        f = open(d["FileName"], "w+")

        f.write("Analytic formula for the resulting magnetic field:\r\n%s\r\n\r\n" % str(B))
        f.write("Data concerning the trajectory of the particle:\r\n")
        f.write("x coordinates of the trajectory:\r\n")
        for x in data["x"]:
            f.write("%f\r\n" % x)
        f.write("y coordinates of the trajectory:\r\n")
        for y in data["y"]:
            f.write("%f\r\n" % y)
        f.write("z coordinates of the trajectory:\r\n")
        for z in data["z"]:
            f.write("%f\r\n" % z)
        f.write("time particle is at a certain position (in sec):\r\n")
        for t in data["t"]:
            f.write("%f\r\n" % t)
        f.write("|v| at a certain position:\r\n")
        for v in data["|v|"]:
            f.write("%f\r\n" % v)
        f.write("|a| at a certain position:\r\n")
        for a in data["|a|"]:
            f.write("%f\r\n" % a)
        f.write("End of file\r\n")

        f.close()

