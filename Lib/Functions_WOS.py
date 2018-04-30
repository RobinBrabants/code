# Contains the numerical method function for the Walk On Spheres method used to evaluate certain electrical fields at specific points used in the Particle::ParticleMove function


from random import random
from Lib.Objects_3D import *
from operator import itemgetter

from sympy.vector import CoordSys3D, Vector
import math


def potential_WOS(electrodes_WOS, Coordinates):
    # see the documentation on how the Walk On Spheres method works (Walk On Spheres)

    max_dist = 10 ** 4  # can be chosen accordingly to the dimensions of the electrode setup
    # If the calculated minimal distance to an electrode exceeds max_dist, it results in choosing a potential equal to V_inf for that iteration
    # Because the chance of reaching an electrode surface then will be too small
    V_inf = 0               # potential on infinity
    max_iter = 400          # if the number of iterations for reaching an electrode surface becomes bigger than max_iter, this step in the numerical method will be skipped
    space = 10 ** (-2)      # can be chosen accordingly to the dimensions of the electrode setup
    # maximum distance of the point from a certain iteration to the surface of an electrode for which the potential of this surface will be taken for that iteration
    iterations = 2000        # number of times the WOS method should be executed before assigning a potential to a point

    # starting values:
    Potential = 0
    iteration = 0
    values = 0          # keeps track on how many times the WOS method was successfully executed to calculate the average potential in a later stage


    while iteration < iterations:
        iter = 0
        Point = Coordinates[:]

        while iter < max_iter:

            distances = []

            for electrode in electrodes_WOS:
                distances.append(electrode.GetClosestDistanceToPoint(Point))

            index, distance = min(enumerate(distances), key=itemgetter(1))

            if distance >= max_dist:
                # if the minimal distance to an electrode surface exceeds max_dist, an infinity potential will be chosen for that iteration:
                Potential += V_inf

                values += 1

                break

            elif distance <= space:
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

        print('%f %%  completed' % (100*(iteration / iterations)), end='\r')

        iteration += 1

    # Take the average
    Potential = Potential / values

    print("\r\n")

    return Potential


def ElectricalField_WOS(electrodes_WOS, Coordinates):
    # see the documentation on how the Walk On Spheres method works (Walk On Spheres)

    L = CoordSys3D('L')
    pi = math.pi

    max_dist = 10 ** 4  # can be chosen accordingly to the dimensions of the electrode setup
    # If the calculated minimal distance to an electrode exceeds max_dist, it results in choosing an electrical field vector equal to E_inf for that iteration
    # Because the chance of reaching an electrode surface then will be too small
    E_inf = Vector.zero     # electrical field on infinity
    max_iter = 400          # if the number of iterations for reaching an electrode surface becomes bigger than max_iter, this step in the numerical method will be skipped
    space = 10 ** (-2)      # can be chosen accordingly to the dimensions of the electrode setup
    # maximum distance of the point from a certain iteration to the surface of an electrode for which the potential of this surface will be taken for that iteration
    iterations = 2000        # number of times the WOS method should be executed before assigning a potential to a point

    # starting values:
    ElectricalField = Vector.zero
    iteration = 0
    values = 0  # keeps track on how many times the WOS method was successfully executed to calculate the average potential in a later stage

    Point_Circle = [0, 0, 0]

    while iteration < iterations:
        iter = 0
        Point = Coordinates[:]


        while iter < max_iter:

            distances = []

            for electrode in electrodes_WOS:
                distances.append(electrode.GetClosestDistanceToPoint(Point))

            index, distance = min(enumerate(distances), key=itemgetter(1))

            if distance >= max_dist:
                # if the minimal distance to an electrode surface exceeds max_dist, an infinity electrical field vector will be chosen for that iteration:
                ElectricalField += E_inf

                values += 1

                break

            elif distance <= space:
                Potential = electrodes_WOS[index].Potential


                r = (Coordinates[0] - Point_Circle[0])*L.i + (Coordinates[1] - Point_Circle[1])*L.j + (Coordinates[2] - Point_Circle[2])*L.k
                er = r.normalize()

                ElectricalField += (3/(r.magnitude()))*Potential*er

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

                if iter == 0:
                    Point_Circle[0] = Point[0]
                    Point_Circle[1] = Point[1]
                    Point_Circle[2] = Point[2]

            iter += 1



        print('%f %%  completed' % (100*(iteration / iterations)), end='\r')

        iteration += 1

    # Take the average
    ElectricalField = ElectricalField / values

    print("\r\n")

    return ElectricalField