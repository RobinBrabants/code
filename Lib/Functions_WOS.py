from scipy.integrate import quad
import xml.etree.ElementTree as ET
from Lib.MagneticFieldComponents import *
from Lib.ElectricFieldComponents import *
from sympy import sqrt, diff, solve, nsolve
from scipy.optimize import fsolve, fmin
import numpy as np
import math
from random import random
from operator import itemgetter
import multiprocessing as mp
from Lib.Elements import *
from multiprocessing import pool
import time


def Get_Data_EWOS(d):
    # function to get the data (list of electrode classes) for the main function
    electrodes = []

    i = 1
    while "Msphere2_" + str(i) in d:
        electrodes.append(Sphere("sphere_" + str(i) , d["Msphere2_" + str(i)], d["Rsphere2_" + str(i)], d["Vsphere2_" + str(i)]))
        i += 1

    i = 1
    while "Pcylinder_" + str(i) in d:
        electrodes.append(Cylinder(d["Pcylinder_" + str(i)], d["Qcylinder_" + str(i)], d["Rcylinder_" + str(i)],
                                 d["Vcylinder_" + str(i)]))
        i += 1

    i = 1
    while "Mcircdisc_" + str(i) in d:
        electrodes.append(CircularDisk(d["Mcircdisc_" + str(i)], d["Rcircdisc_" + str(i)], d["Phicircdisc_" + str(i)],
                                     d["Thetacircdisc_" + str(i)], d["Vcircdisc_" + str(i)]))
        i += 1

    i = 1
    while "Mcircdisc4h_" + str(i) in d:
        electrodes.append(
            CircularDisk4Holes(d["Mcircdisc4h_" + str(i)], d["Mcircdisc4h1_" + str(i)], d["Mcircdisc4h2_" + str(i)],
                               d["Mcircdisc4h3_" + str(i)], d["Mcircdisc4h4_" + str(i)], d["Rcircdisc4h_" + str(i)],
                               d["Rcircdisc4h1_" + str(i)], d["Rcircdisc4h2_" + str(i)],
                               d["Rcircdisc4h3_" + str(i)], d["Rcircdisc4h4_" + str(i)], d["Ncircdisc4h_" + str(i)],
                               d["Vcircdisc4h_" + str(i)]))
        i += 1


    return electrodes


def WalkOnSpheres_potential_slice(d, x, y, z):
    global V_inf, max_dist, maxsteps, space, iterations
    # constanten voor de potential_EWOS functie
    V_inf = 0       #potentiaal op oneindig
    max_dist = 10 ** 8  # arbitrair in te stellen afhankelijk van de dimensies van het probleem. Over deze afstand tot de oppervlakken gaan resulteert in een V_inf-potentiaal (V op oneindig wordt hier immers V_inf gesteld). Kans is immers te klein om nog een oppervlak te bereiken --> inefficiënt.
    maxsteps = 400  # indien er over dit aantal iteraties wordt gegaan, gooien we de iteratie weg.
    space = 10 ** (-2)  # arbitrair in te stellen afhankelijk van de dimensies van het probleem. Afstand tot oppervlak waarop iteratie de potentiaal van het oppervlak aanneemt.
    iterations = 200  # niet te verwarren met maxsteps


    electrodes = Get_Data_EWOS(d)

    Vgrid = np.zeros(shape=(len(z), len(y)))

    factor = len(y)
    number_iterations = len(y) * len(z)

    #multi-processing of the different points in Vgrid (determines automatically how many parallel processes it can run depending on the amount of CPU's in your computer)
    print("WOS process: start")
    pool = mp.Pool()
    results = [pool.apply_async(MultiProcess_WOS, args=([x, j, i], z[::-1].index(i), y.index(j), electrodes, factor, number_iterations)) for i in z[::-1] for j in y]
    pool.close
    pool.join

    list = [p.get() for p in results]

    #list has to be sorted cause processes don't finish in the correct order
    list.sort()
    list = [r[1] for r in list]

    #assignment to the Vgrid
    for k in range(0,len(z)):
        for l in range(0, len(y)):
            Vgrid[k, l] = list[l+factor*k]

    print("WOS process: completed")

    np.set_printoptions(threshold=np.nan, linewidth=500)

    print(Vgrid)

    #write to the file
    f = open("Vgrid", "w+")

    f.write("x\r\n")
    f.write(str(x))
    f.write("\r\n")
    f.write("y\r\n")
    f.write(str(y))
    f.write("\r\n")
    f.write("z\r\n")
    f.write(str(z))
    f.write("\r\n")
    f.write("Vgrid\r\n")
    f.write(str(Vgrid))

    return Vgrid


def MultiProcess_WOS(P, k, l, electrodes, factor, number_iteration):
    potential = potential_EWOS(electrodes, P)

    # to keep track of the process
    print("WOS process:", l + factor * k + 1, "out of", number_iteration, "(",100 * (l + factor * k + 1) / number_iteration, "percent)")

    return (l + factor * k, potential)


def potential_EWOS(electrodes, P):

    V = 0
    k = 0
    values = 0
    pointinelectrode = 0


    for electrode in electrodes:
        if electrode.ispointinoronobject(P):
            V = electrode.potential
            pointinelectrode = 1

    if pointinelectrode == 0:
        while k < iterations:
            j = 0
            P_var = P[:]
            while j < maxsteps:

                dist = []

                for electrode in electrodes:
                    dist.append(electrode.getclosestdistancetopoint(P_var))

                #bol heeft als straal de kleinste afstand

                index, radius = min(enumerate(dist), key=itemgetter(1))
                #als de straal te groot is wordt de ptentiaal gelijk gesteld aan V_inf
                if radius >= max_dist:
                    V += V_inf

                    values += 1

                    break

                elif radius <= space:
                    V += electrodes[index].potential

                    values += 1

                    break

                else:
                    #random vector aanmaken en P_var updaten
                    x = random()-0.5
                    y = random()-0.5
                    z = random()-0.5

                    length = sqrt(x**2+y**2+z**2)
                    vectorx = x/length
                    vectory = y / length
                    vectorz = z / length

                    P_var[0] += vectorx*radius
                    P_var[1] += vectory*radius
                    P_var[2] += vectorz*radius

                j += 1
            k += 1

        #gemiddelde
        V = V / values

    return V


"""def WalkOnSpheres(EWOS, P):
    x_1 = P[0]
    y_1 = P[1]
    z_1 = P[2]

    L = CoordSys3D('L')
    pi = math.pi
    functions = []
    potentials = []
    i = 0

    print(EWOS)

    if EWOS.find("EvalWalkOnSpheres") != -1:
        while EWOS.find("EvalWalkOnSpheres") != -1:
            end = EWOS.find(")")+1

            x, y = sy.symbols('x y')

            functions.append(EWOS[EWOS.find("function")+11:EWOS.find(", potential")])
            potentials.append(EWOS[EWOS.find("potential")+12:EWOS.find(")")])

            if functions[i].find("sphere") != -1:
                functions[i] = functions[i]
                potentials[i] = eval(potentials[i])
            else:
                functions[i] = eval(functions[i])
                potentials[i] = eval(potentials[i])

            EWOS = EWOS[end:]

            i+=1


    # enkel potentiaal berekenen:

    max_i = i

    V_inf = 0
    max_dist = 10**8 #arbitrair in te stellen afhankelijk van de dimensies van het probleem. Over deze afstand tot de oppervlakken gaan resulteert in een V_inf-potentiaal (V op oneindig wordt hier immers V_inf gesteld). Kans is immers te klein om nog een oppervlak te bereiken --> inefficiënt.
    maxsteps = 400 #indien er over dit aantal iteraties wordt gegaan, gooien we de iteratie weg.
    space = 10 ** (-4) #arbitrair in te stellen afhankelijk van de dimensies van het probleem. Afstand tot oppervlak waarop iteratie de potentiaal van het oppervlak aanneemt.
    iterations = 50 #niet te verwarren met maxsteps

    V = 0
    k = 0
    values = 0

    while k < iterations:
        j = 0
        P_var = P[:]
        while j < maxsteps:
            dist = []

            for i in range(0, max_i):
                if functions[i].find("sphere") != -1:
                    M = eval(functions[i][functions[i].find("["):functions[i].find("]") + 1])
                    R = eval(functions[i][functions[i].find(";") + 2:])

                    dist.append(dist_sphere(P_var, M, R))

                else:
                    function = eval(functions[i])
                    dist.append(dist_surface(P_var, function))

            index, radius = min(enumerate(dist), key=itemgetter(1))


            if radius >= max_dist:
                V += V_inf

                values += 1

                break

            # als het de eerste keer kleiner is dan space hebben we in feite de particlemove functie al stopgezet
            if radius <= space:
                potential = potentials[index]

                V += potential

                values += 1

                break

            else:
                phi = 2 * pi * random()
                theta = pi * random()
                vector = radius * ConvertAnglesToVector(phi, theta)

                vectorc = vector.components
                UpdateDictionary(vectorc)

                P_var[0] += vectorc[L.i]
                P_var[1] += vectorc[L.j]
                P_var[2] += vectorc[L.k]

            j += 1

        k += 1
        print(k)

    V = V / values

    print(V)

    return V"""


#functie voor E i.p.v. V, werkt nog niet
def WalkOnSpheres(EWOS, P):
    x_1 = P[0]
    y_1 = P[1]
    z_1 = P[2]

    L = CoordSys3D('L')
    pi = math.pi
    functions = []
    potentials = []
    i = 0

    print(EWOS)

    if EWOS.find("EvalWalkOnSpheres") != -1:
        while EWOS.find("EvalWalkOnSpheres") != -1:
            end = EWOS.find(")")+1

            x, y = sy.symbols('x y')

            functions.append(EWOS[EWOS.find("function")+11:EWOS.find(", potential")])
            potentials.append(EWOS[EWOS.find("potential")+12:EWOS.find(")")])

            if functions[i].find("sphere") != -1:
                functions[i] = functions[i]
                potentials[i] = eval(potentials[i])
            else:
                functions[i] = eval(functions[i])
                potentials[i] = eval(potentials[i])

            EWOS = EWOS[end:]

            i+=1

    max_i = i
    maxsteps = 1000
    space = 10**(-4)


    Point = [0, 0, 0]
    E = Vector.zero
    k = 0
    values = 0



    while k <= 300:
        j = 0
        P_var = P[:]
        while j < 100:
            dist = []

            for i in range(0, max_i):
                if functions[i].find("sphere")!= -1:
                    M = eval(functions[i][functions[i].find("["):functions[i].find("]")+1])
                    R = eval(functions[i][functions[i].find(";")+2:])

                    dist.append(dist_sphere(P_var, M, R))

                else:
                    function = eval(functions[i])
                    dist.append(dist_surface(P_var, function))


            index, radius = min(enumerate(dist), key=itemgetter(1))



            #als het de eerste keer kleiner is dan space hebben we in feite de particlemove functie al stopgezet
            if radius <= space:
                potential = potentials[index]


                r = (Point[0] - P[0]) * L.i + (Point[1] - P[1]) * L.j + (Point[2] - P[2]) * L.k
                er = r.normalize()


                E += (-3/(4*pi*r.magnitude()**3))*potential*er

                values += 1

                break

            else:
                phi = 2*pi*random()
                theta = pi*random()
                vector = radius*ConvertAnglesToVector(phi, theta)


                vectorc = vector.components
                UpdateDictionary(vectorc)

                P_var[0] += vectorc[L.i]
                P_var[1] += vectorc[L.j]
                P_var[2] += vectorc[L.k]

                if j == 0:
                    Point[0] = P_var[0]
                    Point[1] = P_var[1]
                    Point[2] = P_var[2]


            j+=1

        k+=1

    E = E/k

    print(E)

    return E












def dist_sphere(P, M, R):
    from math import sqrt

    x_1 = P[0]
    y_1 = P[1]
    z_1 = P[2]
    x_2 = M[0]
    y_2 = M[1]
    z_2 = M[2]

    distance = sqrt((x_1-x_2)**2+(y_1-y_2)**2+(z_1-z_2)**2)-R

    return distance


def dist_cylinder(P, S, Q, R):
    from math import sqrt


    x = P[0]
    y = P[1]
    z = P[2]
    x_1 = S[0]
    y_1 = S[1]
    z_1 = S[2]
    x_2 = Q[0]
    y_2 = Q[1]
    z_2 = Q[2]

    L = CoordSys3D('L')

    v = (x_2-x_1)*L.i + (y_2-y_1)*L.j + (z_2-z_1)*L.k
    e_v = v.normalize()
    e_v_c = e_v.components
    UpdateDictionary(e_v_c)
    a = e_v_c[L.i]
    b = e_v_c[L.j]
    c = e_v_c[L.k]

    if c != 0:
        if (1/c)*(a*x_1+b*y_1+c*z_1-a*x-b*y)<=z<=(1/c)*(a*x_2+b*y_2+c*z_2-a*x-b*y):
            t = (a*x+b*y+c*z-a*x_1-b*y_1-z_1)/(c*(z_2-z_1)+a*(x_2-x_1)+b*(y_2-y_1))
            x_c = x_1 + t*(x_2-x_1)
            y_c = y_1 + t * (y_2 - y_1)
            z_c = z_1 + t * (z_2 - z_1)
            distance = abs(sqrt((x-x_c)**2+(y-y_c)**2+(z-z_c)**2)-R)
        elif z>(1/c)*(a*x_2+b*y_2+c*z_2-a*x-b*y):
            delta = (x-x_2)*L.i + (y-y_2)*L.j + (z-z_2)*L.k
            distance = sqrt((e_v.dot(delta))**2 + (abs((e_v.cross(delta)).magnitude())-R)**2)
        else:
            delta = (x-x_1)*L.i + (y-y_1)*L.j + (z-z_1)*L.k
            distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - R) ** 2)
    else:
        if b != 0:
            if (1/b)*(a*x_1+b*y_1+c*z_1-a*x-c*z)<=y<=(1/b)*(a*x_2+b*y_2+c*z_2-a*x-c*z):
                t = (a * x + b * y + c * z - a * x_1 - b * y_1 - z_1) / (c * (z_2 - z_1) + a * (x_2 - x_1) + b * (y_2 - y_1))
                x_c = x_1 + t * (x_2 - x_1)
                y_c = y_1 + t * (y_2 - y_1)
                z_c = z_1 + t * (z_2 - z_1)
                distance = abs(sqrt((x - x_c) ** 2 + (y - y_c) ** 2 + (z - z_c) ** 2) - R )
            elif y>(1/b)*(a*x_2+b*y_2+c*z_2-a*x-c*z):
                delta = (x - x_2) * L.i + (y - y_2) * L.j + (z - z_2) * L.k
                distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - R) ** 2)
            else:
                delta = (x - x_1) * L.i + (y - y_1) * L.j + (z - z_1) * L.k
                distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - R) ** 2)
        else:
            if (1/a)*(a*x_1+b*y_1+c*z_1-b*y-c*z)<=x<=(1/a)*(a*x_2+b*y_2+c*z_2-b*y-c*z):
                t = (a * x + b * y + c * z - a * x_1 - b * y_1 - z_1) / (c * (z_2 - z_1) + a * (x_2 - x_1) + b * (y_2 - y_1))
                x_c = x_1 + t * (x_2 - x_1)
                y_c = y_1 + t * (y_2 - y_1)
                z_c = z_1 + t * (z_2 - z_1)
                distance = abs(sqrt((x - x_c) ** 2 + (y - y_c) ** 2 + (z - z_c) ** 2) - R )
            elif x>(1/a)*(a*x_2+b*y_2+c*z_2-b*y-c*z):
                delta = (x - x_2) * L.i + (y - y_2) * L.j + (z - z_2) * L.k
                distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - R) ** 2)
            else:
                delta = (x - x_1) * L.i + (y - y_1) * L.j + (z - z_1) * L.k
                distance = sqrt((e_v.dot(delta)) ** 2 + (abs((e_v.cross(delta)).magnitude()) - R) ** 2)

    return distance


def dist_circdisc(P, M, R, Phi, Theta):
    from math import sqrt

    Phi = math.radians(Phi)
    Theta = math.radians(Theta)
    n = ConvertAnglesToVector(Phi, Theta)

    x = P[0]
    y = P[1]
    z = P[2]
    x_1 = M[0]
    y_1 = M[1]
    z_1 = M[2]

    L = CoordSys3D('L')

    """a = (x-x_1)*L.i + (y-y_1)*L.j + (z-z_1)*L.k
    D = eval(str(sin(acos((n.dot(a))/a.magnitude()))*a.magnitude()))"""


    n_c = n.components
    UpdateDictionary(n_c)
    a = n_c[L.i]
    b = n_c[L.j]
    c = n_c[L.k]
    x_2 = a+x_1
    y_2 = b+y_1
    z_2 = c+z_1
    D = eval(str((((x_2-x_1)*L.i + (y_2-y_1)*L.j + (z_2-z_1)*L.k).cross((x_1-x)*L.i + (y_1-y)*L.j + (z_1-z)*L.k)).magnitude()))


    if D<= R:
        distance = abs(-a*x-b*y-c*z+a*x_1+b*y_1+c*z_1)
    else:
        delta = (x - x_1) * L.i + (y - y_1) * L.j + (z - z_1) * L.k
        distance = sqrt((n.dot(delta)) ** 2 + (abs((n.cross(delta)).magnitude()) - R) ** 2)

    return distance


def dist_circdisc_4holes(P, M, M1, M2, M3, M4, R, R1, R2, R3, R4, n_vec):
    from math import sqrt

    L = CoordSys3D('L')
    n = n_vec[0]*L.i + n_vec[1]*L.j + n_vec[2]*L.k

    x = P[0]
    y = P[1]
    z = P[2]
    x_1 = M[0]
    y_1 = M[1]
    z_1 = M[2]
    #holes:
    x1 = M1[0]
    y1 = M1[1]
    z1 = M1[2]
    x2 = M2[0]
    y2 = M2[1]
    z2 = M2[2]
    x3 = M3[0]
    y3 = M3[1]
    z3 = M3[2]
    x4 = M4[0]
    y4 = M4[1]
    z4 = M4[2]

    L = CoordSys3D('L')

    a = n_vec[0]
    b = n_vec[1]
    c = n_vec[2]
    x_2 = a+x_1
    y_2 = b+y_1
    z_2 = c+z_1
    D = eval(str((((x_2-x_1)*L.i + (y_2-y_1)*L.j + (z_2-z_1)*L.k).cross((x_1-x)*L.i + (y_1-y)*L.j + (z_1-z)*L.k)).magnitude()))


    if D<= R:
        x_2 = a + x1
        y_2 = b + y1
        z_2 = c + z1
        D1 = eval(str((((x_2-x1)*L.i + (y_2-y1)*L.j + (z_2-z1)*L.k).cross((x1-x)*L.i + (y1-y)*L.j + (z1-z)*L.k)).magnitude()))
        x_2 = a + x2
        y_2 = b + y2
        z_2 = c + z2
        D2 = eval(str((((x_2-x2)*L.i + (y_2-y2)*L.j + (z_2-z2)*L.k).cross((x2-x)*L.i + (y2-y)*L.j + (z2-z)*L.k)).magnitude()))
        x_2 = a + x3
        y_2 = b + y3
        z_2 = c + z3
        D3 = eval(str((((x_2-x3)*L.i + (y_2-y3)*L.j + (z_2-z3)*L.k).cross((x3-x)*L.i + (y3-y)*L.j + (z3-z)*L.k)).magnitude()))
        x_2 = a + x4
        y_2 = b + y4
        z_2 = c + z4
        D4 = eval(str((((x_2-x4)*L.i + (y_2-y4)*L.j + (z_2-z4)*L.k).cross((x4-x)*L.i + (y4-y)*L.j + (z4-z)*L.k)).magnitude()))

        if D1 <= R1:
            delta = (x - x1) * L.i + (y - y1) * L.j + (z - z1) * L.k
            distance = sqrt((n.dot(delta)) ** 2 + (abs((n.cross(delta)).magnitude()) - R1) ** 2)
        elif D2 <= R2:
            delta = (x - x2) * L.i + (y - y2) * L.j + (z - z2) * L.k
            distance = sqrt((n.dot(delta)) ** 2 + (abs((n.cross(delta)).magnitude()) - R2) ** 2)
        elif D3 <= R3:
            delta = (x - x3) * L.i + (y - y3) * L.j + (z - z3) * L.k
            distance = sqrt((n.dot(delta)) ** 2 + (abs((n.cross(delta)).magnitude()) - R3) ** 2)
        elif D4 <= R4:
            delta = (x - x4) * L.i + (y - y4) * L.j + (z - z4) * L.k
            distance = sqrt((n.dot(delta)) ** 2 + (abs((n.cross(delta)).magnitude()) - R4) ** 2)
        else:
            distance = abs(-a*x-b*y-c*z+a*x_1+b*y_1+c*z_1)
    else:
        delta = (x - x_1) * L.i + (y - y_1) * L.j + (z - z_1) * L.k
        distance = sqrt((n.dot(delta)) ** 2 + (abs((n.cross(delta)).magnitude()) - R) ** 2)

    return distance


def dist_surface(P, surface):
    x_1 = P[0]
    y_1 = P[1]
    z_1 = P[2]

    x, y = sy.symbols('x, y')

    def f(t):
        x, y = t
        return eval(str(surface))

    x, y = fmin(f, [x_1, y_1], disp=False)

    z = eval(str(surface))

    radius = sqrt((x_1-x)**2+(y_1-y)**2+(z_1-z)**2)

    return radius
