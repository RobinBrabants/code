# Contains the numerical method function for the Walk On Spheres method used to evaluate certain electrical fields at specific points used in the Particle::ParticleMove function


from random import random
from Lib.Objects_3D import *
from operator import itemgetter

from sympy.vector import CoordSys3D, Vector
import math


def potential_WOS(electrodes_WOS, Coordinates, d):
    # see the documentation on how the Walk On Spheres method works (Walk On Spheres)


    V_inf = 0               # potential on infinity

    # starting values:
    Potential = 0
    iteration = 0
    values = 0          # keeps track on how many times the WOS method was successfully executed to calculate the average potential in a later stage


    while iteration < d["Iterations"]:
        iter = 0
        Point = Coordinates[:]

        while iter < d["MaximumIterations"]:

            distances = []

            for electrode in electrodes_WOS:
                distances.append(electrode.DistanceToObject(Point))

            index, distance = min(enumerate(distances), key=itemgetter(1))

            if distance >= d["MaximumDistance"]:
                # if the minimal distance to an electrode surface exceeds d["MaximumDistance"], an infinity potential will be chosen for that iteration:
                Potential += V_inf

                values += 1

                break

            elif distance <= d["Gap"]:
                Potential += electrodes_WOS[index].Potential

                values += 1

                break

            else:
                # make random vector and update the next point of the iteration:
                x = random()-0.5
                y = random()-0.5
                z = random()-0.5

                length = sqrt(x**2+y**2+z**2)
                vector_x = x / length
                vector_y = y / length
                vector_z = z / length

                Point[0] += vector_x*distance
                Point[1] += vector_y*distance
                Point[2] += vector_z*distance

            iter += 1

        print('%f %%  completed' % (100*(iteration / d["Iterations"])), end='\r')

        iteration += 1

    # Take the average
    Potential = Potential / values

    print("\r\n")

    return Potential


def ElectricalField_WOS(electrodes_WOS, Coordinates, d):
    # see the documentation on how the Walk On Spheres method works (Walk On Spheres)

    L = CoordSys3D('L')

    E_inf = Vector.zero             # electrical field on infinity

    # starting values:
    ElectricalField = Vector.zero
    iteration = 0
    values = 0                      # keeps track on how many times the WOS method was successfully executed to calculate the average electrical field in a later stage
    Point_Circle = [0, 0, 0]

    while iteration < d["Iterations"]:
        iter = 0
        Point = Coordinates[:]


        while iter < d["MaximumIterations"]:

            distances = []

            for electrode in electrodes_WOS:
                distances.append(electrode.DistanceToObject(Point))

            index, distance = min(enumerate(distances), key=itemgetter(1))

            if distance >= d["MaximumDistance"]:
                # if the minimal distance to an electrode surface exceeds d["MaximumDistance"], an infinity electrical field vector will be chosen for that iteration:
                ElectricalField += E_inf

                values += 1

                break

            elif distance <= d["Gap"]:
                Potential = electrodes_WOS[index].Potential

                # test WOS:
                r1 = (Coordinates[0] - Point_Circle[0])*L.i + (Coordinates[1] - Point_Circle[1])*L.j + (Coordinates[2] - Point_Circle[2])*L.k
                r2 = (Coordinates[0] - Point[0]) * L.i + (Coordinates[1] - Point[1]) * L.j + (Coordinates[2] - Point[2]) * L.k

                # Point_Cirlce or Point itself dunno

                er = r2.normalize()

                ElectricalField += (1/(r2.magnitude()))*Potential*er
                #ElectricalField += Potential * er

                values += 1

                break

            else:
                # make random vector and update the next point of the iteration:
                u = random()
                v = random()

                theta = 2* math.pi * u
                phi = math.acos(2*v-1)


                Point[0] += distance*math.cos(theta)*math.sin(phi)
                Point[1] += distance*math.sin(theta)*math.sin(phi)
                Point[2] += distance*math.cos(phi)

                if iter == 0:
                    Point_Circle[0] = Point[0]
                    Point_Circle[1] = Point[1]
                    Point_Circle[2] = Point[2]

            iter += 1



        #print('%f %%  completed' % (100*(iteration / d["Iterations"])), end='\r')

        iteration += 1

    # Take the average
    ElectricalField = ElectricalField / values

    #print("\r\n")

    return ElectricalField