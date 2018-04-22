# Particle class with a few functions (see class itself)


from Lib.Functions import *
from sympy.vector import CoordSys3D, Vector
import sympy as sy
import matplotlib.pyplot as plt


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
        self.Theta1 = Theta1                    # measured from the x-axis in anti-clockwise direction until the projection of the velocity vector on the xy-plane (0-360 degree)
        self.Phi1 = Phi1                        # measured from the projection of the velocity vector on the xy-plane in anti-clockwise direction until the velocity vector (0-360 degree)
        self.Acceleration = Acceleration        # acceleration of the particle which is the result of an external force and will all the time be added to the acceleration caused by the electric and magnetic fields (in m/s²)
        self.Theta2 = Theta2                    # measured from the x-axis in anti-clockwise direction until the projection of the acceleration vector on the xy-plane (0-360 degree)
        self.Phi2 = Phi2                        # measured from the projection of the velocity vector on the xy-plane in anti-clockwise direction until the acceleration vector (0-360 degree)
        self.Position = Position                # starting position particle

    def ParticleMove(self, B_analytic, E_analytic, electrodes, electrodes_WOS, d):
        # calculates the trajectory of a particle in predetermined electric and magnetic fields
        # see the documentation on how the trajectory of the particle is determined (Beweging deeltje in magneetveld)

        print("The trajectory of the particle is now being calculated")

        L = CoordSys3D('L')
        x, y, z = sy.symbols('x y z')

        Coordinates = self.Position
        Speed = GetVector(self.Velocity, self.Theta1,  self.Phi1)
        Acceleration_constant = GetVector(self.Acceleration, self.Theta2, self.Phi2)
        E, B = GetFields(Coordinates, Speed, B_analytic, E_analytic, electrodes_WOS)                    # Step 1, 2
        Acceleration = (self.Charge / self.Mass) * ((Speed.cross(B)) + E) + Acceleration_constant       # Step 3
        time = 0

        Trajectory = {"t": ["time"], "x": ["x coord"], "y": ["y coord"], "z": ["z coord"], "v": ["velocity"], "a": ["acceleration"], "E": ["electric field"], "B": [" magneticfield"]}
        # Dictionary in which all the data about the calculated trajectory will be held
        Trajectory["t"].append(time)
        Trajectory["x"].append(Coordinates[0])
        Trajectory["y"].append(Coordinates[1])
        Trajectory["z"].append(Coordinates[2])
        Trajectory["v"].append(Speed)
        Trajectory["a"].append(Acceleration)
        Trajectory["E"].append(E)
        Trajectory["B"].append(B)

        while time <= d["timelimit"]:           # maximum execute time cannot be exceeded
            # the testpoint is the first estimate of the next point which will later be corrected

            # Step 4
            vComponents = Speed.components
            UpdateDictionary(vComponents)
            Coordinates_testpoint = [Coordinates[0] + vComponents[L.i] * d["timesteps"], Coordinates[1] + vComponents[L.j] * d["timesteps"], Coordinates[2] + vComponents[L.k] * d["timesteps"]]

            # Step 5
            aComponents = Acceleration.components
            UpdateDictionary(aComponents)
            Speed_testpoint = (vComponents[L.i] + aComponents[L.i] * d["timesteps"]) * L.i + (vComponents[L.j] + aComponents[L.j] * d["timesteps"]) * L.j + (vComponents[L.k] + aComponents[L.k] * d["timesteps"]) * L.k
            # best estimate for the speed at the testpoint available
            E_testpoint, B_testpoint = GetFields(Coordinates_testpoint, Speed_testpoint, B_analytic, E_analytic, electrodes_WOS)
            Acceleration_testpoint = (self.Charge / self.Mass) * ((Speed_testpoint.cross(B_testpoint)) + E_testpoint) + Acceleration_constant

            # Step 6
            Acceleration_testpoint_average = (Acceleration + Acceleration_testpoint)/2

            # Step 7
            aComponents_testpoint_average = Acceleration_testpoint_average.components
            UpdateDictionary(aComponents_testpoint_average)
            Speed_testpoint = (vComponents[L.i] + aComponents_testpoint_average[L.i] * d["timesteps"]) * L.i + (vComponents[L.j] + aComponents_testpoint_average[L.j] * d["timesteps"]) * L.j + (vComponents[L.k] + aComponents_testpoint_average[L.k] * d["timesteps"]) * L.k
            # better estimate for the speed at the testpoint

            # Step 8
            Speed_testpoint_average = (Speed + Speed_testpoint)/2

            # Step 9
            vComponents_testpoint_average = Speed_testpoint_average.components
            UpdateDictionary(vComponents_testpoint_average)
            Coordinates = [Coordinates[0] + vComponents_testpoint_average[L.i] * d["timesteps"], Coordinates[1] + vComponents_testpoint_average[L.j] * d["timesteps"], Coordinates[2] + vComponents_testpoint_average[L.k] * d["timesteps"]]

            # Step 10, 11
            E, B = GetFields(Coordinates, Speed_testpoint, B_analytic, E_analytic, electrodes_WOS)
            # Speed_testpoint is a decent estimate for the speed at the next point

            # Step 12
            Acceleration_nextpoint = (self.Charge / self.Mass) * ((Speed_testpoint.cross(B)) + E) + Acceleration_constant
            # Speed_testpoint is a decent estimate for the speed at the next point

            # Step 13
            Acceleration_average = (Acceleration + Acceleration_nextpoint) / 2

            # Step 14
            aComponents_average = Acceleration_average.components
            UpdateDictionary(aComponents_average)
            Speed = (vComponents[L.i] + aComponents_average[L.i] * d["timesteps"]) * L.i + (vComponents[L.j] + aComponents_average[L.j] * d["timesteps"]) * L.j + (vComponents[L.k] + aComponents_average[L.k] * d["timesteps"]) * L.k

            Acceleration = Acceleration_nextpoint

            time += d["t"]



            if Coordinates[0] < d["xmin"] or Coordinates[0] > d["xmax"] or Coordinates[1] < d["ymin"] or Coordinates[1] > d["ymax"] or Coordinates[2] < d["zmin"] or Coordinates[2] > d["zmax"]:
                print("The particle went out of the predetermined box")
                break


            for electrode in electrodes, electrodes_WOS:
                if electrode.IsPointInObject() != 0:
                    print("The particle collided with: " + electrode.name)
                    Trajectory["collision"] = electrode
                    break


            Trajectory["t"].append(time)
            Trajectory["x"].append(Coordinates[0])
            Trajectory["y"].append(Coordinates[1])
            Trajectory["z"].append(Coordinates[2])
            Trajectory["v"].append(Speed)
            Trajectory["a"].append(Acceleration)
            Trajectory["E"].append(E)
            Trajectory["B"].append(B)


        print("The calculations for the trajectory of the particle have been successfully finished")

        return Trajectory

    def WriteFile(self, E_analytic, B_analytic, Trajectory, d):
        # function which writes all the data concerning the calculated trajectory to a file

        print("Writing data from the calculated trajectory of the particle to: " + d["FileName"])

        f = open(d["FileName"], "wb")

        f.write("This file includes the data from the calculated trajectory of the particle: %s, %s" % (self.name, self.Type))

        f.write("Analytic formula for the resulting electric field:\r\n%s\r\n\r\n" % str(E_analytic))
        f.write("Analytic formula for the resulting magnetic field:\r\n%s\r\n\r\n" % str(B_analytic))
        f.write("REMARK: This does not include the electric field produced by the WOS objects \r\n\r\n")

        f.write("Data concerning the trajectory of the particle:\r\n")

        for val in zip( Trajectory["t"], Trajectory["x"], Trajectory["y"], Trajectory["z"], Trajectory["v"], Trajectory["a"], Trajectory["E"], Trajectory["B"]):
            f.write('{}, {}, {}, {}, {}, {}, {}, {}\r\n'.format(val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]))

        if "collision" in Trajectory:
            electrode = Trajectory["collision"]
            f.write("Particle collided with: %s\r\n\r\n" % str(electrode.name))         # says if there was a collision

        f.write("End of file\r\n")

        f.close()

        print(d["FileName"] + "written")

    def PlotTrajectory(self, Trajectory):
        # funtion which plots the trajectory of the particle

        print("Plotting the trajectory of the particle")

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        x = Trajectory["x"]
        y = Trajectory["y"]
        z = Trajectory["z"]
        t = Trajectory["t"]

        xs = np.array(x)
        ys = np.array(y)
        zs = np.array(z)
        ts = np.array(t)

        ax.plot(x, y, z)
        #ax.scatter maybe also good

        len_t = len(t) - 1
        ind_pos = []

        for i in range(0, 11):
            ind_pos.append(0 + i * (int(len_t / 10)))           #time labels

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
        fig.suptitle('trajectory of particle: ' + self.name + ', ' + self.Type, fontsize=14, fontweight='bold')
        ax.set_title('time in seconds', style='italic', fontsize=8, bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 7})

        plt.show(block="False")



