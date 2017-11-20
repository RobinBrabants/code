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
from multiprocessing import pool


def ConvertAnglesToVector(Phi, Theta):
    print("test")
    pi = math.pi
    a = math.tan(Theta)/math.sqrt(1 + math.tan(Theta)**2)

    if 0 <= Phi <= pi / 2 or (3 / 2) * pi <= Phi < 2 * pi:
        ex = abs(math.cos(Phi) * a)
    else:
        ex = -abs(math.cos(Phi) * a)

    if 0<=Phi<=pi:
        ey = abs(math.sqrt(a**2-ex**2))
    else:
        ey = -abs(math.sqrt(a ** 2 - ex ** 2))

    if math.tan(Theta) != 0:
        if 0<=Theta<=pi/2:
            ez = abs(a/math.tan(Theta))
        else:
            ez = -abs(a / math.tan(Theta))
    else:
        ez = 1

    L = CoordSysCartesian('L')

    return ex*L.i + ey*L.j + ez*L.k


def ConvertVectorToAngles(v):

    vcomponents = v.components
    UpdateDictionary(vcomponents)

    L = CoordSysCartesian('L')

    ex = vcomponents[L.i]
    ey = vcomponents[L.j]
    ez = vcomponents[L.k]

    a = math.sqrt(ey**2 + ex**2)
    Phi = math.degrees(math.acos(ex/a))
    Theta = math.degrees(math.atan(a/ez))

    return Phi, Theta


def UpdateDictionary(Dict):
    L = CoordSysCartesian('L')
    if L.i not in Dict:
        Dict.setdefault(L.i, 0)
    if L.j not in Dict:
        Dict.setdefault(L.j, 0)
    if L.k not in Dict:
        Dict.setdefault(L.k, 0)


def EvaluateField(B, S):

    x, y, z, Phi2, t = sy.symbols('x y z Phi2 t')
    L = CoordSysCartesian('L')


    Bcomponents = B.components
    UpdateDictionary(Bcomponents)

    for basis in [L.i, L.j, L.k]:

        sum_integrals = 0


        if str(Bcomponents[basis]).find("Integral(") != -1:
            while str(str(Bcomponents[basis])).find("Integral(") != -1:
                str(B)
                """This solved it somehow B stays the same now throughout the process"""
                end1 = str(str(Bcomponents[basis])).find("Integral(")

                integrand = str(Bcomponents[basis])[str(Bcomponents[basis]).find("Integral(") + 9:str(Bcomponents[basis]).find(", ")]
                integrand = sy.sympify(integrand)
                integrand = integrand.subs([(x, S[0]), (y, S[1]), (z, S[2])])

                info1 = str(Bcomponents[basis])[str(Bcomponents[basis]).find(",") + 1:]
                dx = info1[info1.find("(") + 1:info1.find(",")]
                info2 = info1[info1.find(",") + 2:]
                a = info2[:info2.find(",")]
                b = info2[info2.find(",") + 2:info2.find("))")]

                begin2 = info2.find("))")


                if dx == "Phi2":
                    sum_integrals += quad(lambda Phi2: eval(str(integrand)), eval(a), eval(b))[0]
                if dx == "t":
                    sum_integrals += sy.integrate(eval(str(integrand)), (t,  eval(a), eval(b)))

                Bcomponents[basis] = str(Bcomponents[basis])[:end1] + info2[begin2 + 3:]

            if sum_integrals != 0:
                Bcomponents[basis] += str(sum_integrals)

            Bcomponents[basis] = eval(Bcomponents[basis])


    Bdic = Bcomponents[L.i]*L.i + Bcomponents[L.j]*L.j + Bcomponents[L.k]*L.k

    if Bdic != Vector.zero:
        Beval = Bdic.evalf(subs={x: S[0], y: S[1], z: S[2]})
    else:
        Beval = Vector.zero

    return Beval


def WalkOnSpheres_potential_slice(EWOS, x, y, z):

    functions, potentials = Get_Data_EWOS(EWOS)

    print(functions)

    Vgrid = np.zeros(shape=(len(z), len(y)))

    factor = len(y)
    number_iterations = len(y) * len(z)

    #multi-processing of the different points in Vgrid (determines automatically how many parallel processes it can run depending on the amount of CPU's in your computer)
    pool = mp.Pool()
    results = [pool.apply_async(MultiProcess_WOS, args=([x, j, i], z[::-1].index(i), y.index(j), functions, potentials, factor, number_iterations)) for i in z[::-1] for j in y]
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

    print("WOS: completed")

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


def MultiProcess_WOS(P, k, l, functions, potentials, factor, number_iteration):
    potential = potential_EWOS(functions, potentials, P)

    #to keep track of the process
    print("WOS process:",l + factor * k+1,"out of", number_iteration,"(",100*(l + factor * k+1)/number_iteration,"percent)")

    return (l + factor * k, potential)


def Get_Data_EWOS(EWOS):
    # function to get the data (surfaces and their potentials) for the main function
    functions = []
    potentials = []
    i = 0

    if EWOS.find("EvalWalkOnSpheres") != -1:
        while EWOS.find("EvalWalkOnSpheres") != -1:
            end = EWOS.find(")")+1

            functions.append(EWOS[EWOS.find("function")+11:EWOS.find(", potential")])
            potentials.append(EWOS[EWOS.find("potential")+12:EWOS.find(")")])


            if functions[i].find("sphere") != -1:
                functions[i] = functions[i]
                potentials[i] = eval(potentials[i])
            elif functions[i].find("cylinder") != -1:
                functions[i] = functions[i]
                potentials[i] = eval(potentials[i])
            elif functions[i].find("circulardisk") != -1:
                functions[i] = functions[i]
                potentials[i] = eval(potentials[i])
            elif functions[i].find("circdisk4holes") != -1:
                functions[i] = functions[i]
                potentials[i] = eval(potentials[i])
            else:
                functions[i] = eval(functions[i])
                potentials[i] = eval(potentials[i])

            EWOS = EWOS[end:]

            i+=1


    return functions, potentials


def potential_EWOS(functions,potentials,P):

    max_i = len(functions)

    V_inf = 0
    max_dist = 10 ** 8  # arbitrair in te stellen afhankelijk van de dimensies van het probleem. Over deze afstand tot de oppervlakken gaan resulteert in een V_inf-potentiaal (V op oneindig wordt hier immers V_inf gesteld). Kans is immers te klein om nog een oppervlak te bereiken --> inefficiënt.
    maxsteps = 400  # indien er over dit aantal iteraties wordt gegaan, gooien we de iteratie weg.
    space = 10 ** (-2)  # arbitrair in te stellen afhankelijk van de dimensies van het probleem. Afstand tot oppervlak waarop iteratie de potentiaal van het oppervlak aanneemt.
    iterations = 200  # niet te verwarren met maxsteps

    V = 0
    k = 0
    values = 0

    while k < iterations:
        j = 0
        P_var = P[:]
        while j < maxsteps:
            dist = []
            #bepalen van de afstand vanuit P_var tot aan alle oppervlakken
            for i in range(0, max_i):
                if functions[i].find("sphere") != -1:
                    M = eval(functions[i][functions[i].find("["):functions[i].find("]") + 1])
                    R = eval(functions[i][functions[i].find(";") + 2:])
                    dist.append(dist_sphere(P_var, M, R))


                elif functions[i].find("cylinder") != -1:
                    P = eval(functions[i][functions[i].find("["):functions[i].find("]") + 1])
                    Q = eval(functions[i][functions[i].find(";") + 2:functions[i].find(", radius")])
                    R = eval(functions[i][functions[i].find("radius") + 9:])

                    dist.append(dist_cylinder(P_var, P, Q, R))

                elif functions[i].find("circulardisk") != -1:
                    M = eval(functions[i][functions[i].find("["):functions[i].find("]") + 1])
                    R = eval(functions[i][functions[i].find(";") + 2:functions[i].find(", orientation")])
                    Phi = eval(functions[i][functions[i].find("phi") + 4:functions[i].find(", theta")])
                    Theta = eval(functions[i][functions[i].find("theta") + 6:])

                    dist.append(dist_circdisc(P_var, M, R, Phi, Theta))

                elif functions[i].find("circdisk4holes") != -1:
                    M = eval(functions[i][functions[i].find("["):functions[i].find("]") + 1])
                    R = eval(functions[i][functions[i].find(";") + 2:functions[i].find(", orientation")])
                    N = eval(functions[i][functions[i].find("orientation") + 14:functions[i].find(", holes")])
                    M1 = eval(functions[i][functions[i].find("M1") + 3:functions[i].find(", R1")])
                    R1 = eval(functions[i][functions[i].find("R1") + 3:functions[i].find("; M2")])
                    M2 = eval(functions[i][functions[i].find("M2") + 3:functions[i].find(", R2")])
                    R2 = eval(functions[i][functions[i].find("R2") + 3:functions[i].find("; M3")])
                    M3 = eval(functions[i][functions[i].find("M3") + 3:functions[i].find(", R3")])
                    R3 = eval(functions[i][functions[i].find("R3") + 3:functions[i].find("; M4")])
                    M4 = eval(functions[i][functions[i].find("M4") + 3:functions[i].find(", R4")])
                    R4 = eval(functions[i][functions[i].find("R4") + 3:])

                    dist.append(dist_circdisc_4holes(P_var, M, M1, M2, M3, M4, R, R1, R2, R3, R4, N))

                else:
                    function = eval(functions[i])
                    dist.append(dist_surface(P_var, function))

            #bol heeft als straal de kleinste afstand
            index, radius = min(enumerate(dist), key=itemgetter(1))

            #als de straal te groot is wordt de ptentiaal gelijk gesteld aan V_inf
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
                #random vector aanmaken en P_var updaten
                x = random()-0.5
                y = random() - 0.5
                z = random() - 0.5


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

    L = CoordSysCartesian('L')
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

    L = CoordSysCartesian('L')
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

    L = CoordSysCartesian('L')

    v = (x_2-x_1)*L.i + (y_2-y_1)*L.j + (z_2-z_1)*L.k
    e_v = v.normalize()
    e_v_c = e_v.components
    UpdateDictionary(e_v_c)
    a = e_v_c[L.i]
    b = e_v_c[L.j]
    c = e_v_c[L.k]

    if c != 0:
        if (1/c)*(a*x_1+b*y_1+c*z_1-a*x-b*y)<=z<=(1/c)*(a*x_2+b*y_2+c*z_2-a*x-b*y):
            t = (a*x+b*y+c*z)/(c*(z_2-z_1)+a*(x_2-x_1)+b*(y_2-y_1))
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
                t = (a * x + b * y + c * z) / (c * (z_2 - z_1) + a * (x_2 - x_1) + b * (y_2 - y_1))
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
                t = (a * x + b * y + c * z) / (c * (z_2 - z_1) + a * (x_2 - x_1) + b * (y_2 - y_1))
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

    L = CoordSysCartesian('L')

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

    L = CoordSysCartesian('L')
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

    L = CoordSysCartesian('L')

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


def getvector(value, Theta, Phi):

    eI = ConvertAnglesToVector(math.radians(Phi), math.radians(Theta))

    return eI*value


def ReadXml():
    tree = ET.parse("DataFile.xml")
    root = tree.getroot()

    d = {}
    ParticleDictionary = {"electron": [9.10938356*(10**(-31)), -1.6021766208*(10**(-19))], "proton": [1.672621898*(10**(-27)), 1.6021766208*10**(-19)]}


    i = 1
    for StraightConductor in root.iter("StraightConductor"):
        if StraightConductor.attrib["status"] == "enabled":
            CoordinatesPoint1 = StraightConductor.find("CoordinatesPoint1")
            d["P_" + str(i)] = [eval(CoordinatesPoint1.find("x").text), eval(CoordinatesPoint1.find("y").text), eval(CoordinatesPoint1.find("z").text)]
            CoordinatesPoint2 = StraightConductor.find("CoordinatesPoint2")
            d["Q_" + str(i)] = [eval(CoordinatesPoint2.find("x").text), eval(CoordinatesPoint2.find("y").text), eval(CoordinatesPoint2.find("z").text)]
            d["I_" + str(i)] = eval(StraightConductor.find("Current").text)
            i += 1

    i = 1
    for StraightConductorCollection in root.iter("StraightConductorCollection"):
        if StraightConductorCollection.attrib["status"] == "enabled":
            amount = StraightConductorCollection.attrib["amount"]
            for j in range(1, int(amount) + 1):
                CoordinatesPoint = StraightConductorCollection.find("CoordinatesPoint" +str(j))
                d[str(i) + "S_" + str(j)] = [eval(CoordinatesPoint.find("x").text), eval(CoordinatesPoint.find("y").text),
                                eval(CoordinatesPoint.find("z").text)]
            d["Icollection_" + str(i)] = eval(StraightConductorCollection.find("Current").text)
            i += 1

    i = 1
    for RectangularCoil in root.iter("RectangularCoil"):
        if RectangularCoil.attrib["status"] == "enabled":
            StartingPoint = RectangularCoil.find("StartingPoint")
            d["w_" + str(i)] = eval(RectangularCoil.find("Width").text)
            d["l_" + str(i)] = eval(RectangularCoil.find("Length").text)
            d["h_" + str(i)] = eval(RectangularCoil.find("Heigth").text)
            d["N_" + str(i)] = eval(RectangularCoil.find("Windings").text)
            d["SP_" + str(i)] = [eval(StartingPoint.find("x").text), eval(StartingPoint.find("y").text),
                                 eval(StartingPoint.find("z").text)]
            d["Phireccoil_" + str(i)] = eval(RectangularCoil.find("Phi").text)
            d["Thetareccoil_" + str(i)] = eval(RectangularCoil.find("Theta").text)
            d["Psireccoil_" + str(i)] = eval(RectangularCoil.find("Psi").text)
            d["Ireccoil_" + str(i)] = eval(RectangularCoil.find("Current").text)
            d["begin_" + str(i)] = eval(RectangularCoil.find("begin").text)
            i += 1

    i = 1
    for CircularConductor in root.iter("CircularConductor"):
        if CircularConductor.attrib["status"] == "enabled":
            CoordinatesCentre = CircularConductor.find("CoordinatesCentre")
            d["M_" + str(i)] = [eval(CoordinatesCentre.find("x").text), eval(CoordinatesCentre.find("y").text), eval(CoordinatesCentre.find("z").text)]
            d["R_" + str(i)] = eval(CircularConductor.find("Radius").text)
            Orientation = CircularConductor.find("Orientation")
            d["Phi_" + str(i)] = math.radians(eval(Orientation.find("Phi").text))
            d["Theta_" + str(i)] = math.radians(eval(Orientation.find("Theta").text))
            d["Icircle_" + str(i)] = eval(CircularConductor.find("Current").text)
            i += 1

    i = 1

    for BentConductor in root.iter("BentConductor"):
        if BentConductor.attrib["status"] == "enabled":
            CoordinatesCentre = BentConductor.find("CoordinatesCentre")
            Interval = BentConductor.find("Interval")
            d["Mbent_" + str(i)] = [eval(CoordinatesCentre.find("x").text), eval(CoordinatesCentre.find("y").text), eval(CoordinatesCentre.find("z").text)]
            d["Rbent_" + str(i)] = eval(BentConductor.find("Radius").text)
            Orientation = BentConductor.find("Orientation")
            d["Phibent_" + str(i)] = math.radians(eval(Orientation.find("Phi").text))
            d["Thetabent_" + str(i)] = math.radians(eval(Orientation.find("Theta").text))
            d["Ibent_" + str(i)] = eval(BentConductor.find("Current").text)
            d["Interval_" + str(i)] = [eval(Interval.find("from").text), eval(Interval.find("to").text)]
            i += 1

    i = 1

    for CircularCoil in root.iter("CircularCoil"):
        if CircularCoil.attrib["status"] == "enabled":
            CoordinatesCentreCoil = CircularCoil.find("CoordinatesCentre")
            d["hcoil_" + str(i)] = eval(CircularCoil.find("Heigth").text)
            d["Ncoil_" + str(i)] = eval(CircularCoil.find("Windings").text)
            d["Mcoil_" + str(i)] = [eval(CoordinatesCentreCoil.find("x").text), eval(CoordinatesCentreCoil.find("y").text),
                                 eval(CoordinatesCentreCoil.find("z").text)]
            d["Rcoil_" + str(i)] = eval(CircularCoil.find("Radius").text)
            d["Phicoil_" + str(i)] = eval(CircularCoil.find("Phi").text)
            d["Thetacoil_" + str(i)] = eval(CircularCoil.find("Theta").text)
            d["Icoil_" + str(i)] = eval(CircularCoil.find("Current").text)
            d["begincoil_" + str(i)] = eval(CircularCoil.find("begin").text)
            i += 1

    i = 1
    for Sphere in root.iter("Sphere"):
        if Sphere.attrib["status"] == "enabled":
            CoordinatesCentre = Sphere.find("CoordinatesCentre")
            d["Msphere_" + str(i)] = [eval(CoordinatesCentre.find("x").text), eval(CoordinatesCentre.find("y").text), eval(CoordinatesCentre.find("z").text)]
            d["Rsphere_" + str(i)] = eval(Sphere.find("Radius").text)
            d["Qsphere_" + str(i)] = eval(Sphere.find("Charge").text)
            i += 1

    i = 1
    for Sphere2 in root.iter("Sphere2"):
        if Sphere2.attrib["status"] == "enabled":
            CoordinatesCentre = Sphere2.find("CoordinatesCentre")
            d["Msphere2_" + str(i)] = [eval(CoordinatesCentre.find("x").text), eval(CoordinatesCentre.find("y").text), eval(CoordinatesCentre.find("z").text)]
            d["Rsphere2_" + str(i)] = eval(Sphere2.find("Radius").text)
            d["Vsphere2_" + str(i)] = eval(Sphere2.find("Potential").text)

            i += 1
    i = 1
    for Cylinder in root.iter("Cylinder"):
        if Cylinder.attrib["status"] == "enabled":
            CoordinatesPoint1 = Cylinder.find("CoordinatesPoint1")
            d["Pcylinder_" + str(i)] = [eval(CoordinatesPoint1.find("x").text), eval(CoordinatesPoint1.find("y").text), eval(CoordinatesPoint1.find("z").text)]
            CoordinatesPoint2 = Cylinder.find("CoordinatesPoint2")
            d["Qcylinder_" + str(i)] = [eval(CoordinatesPoint2.find("x").text), eval(CoordinatesPoint2.find("y").text), eval(CoordinatesPoint2.find("z").text)]
            d["Rcylinder_" + str(i)] = eval(Cylinder.find("Radius").text)
            d["Vcylinder_" + str(i)] = eval(Cylinder.find("Potential").text)
            i += 1

    i = 1
    for CircularDisc in root.iter("CircularDisc"):
        if CircularDisc.attrib["status"] == "enabled":
            CoordinatesCentre = CircularDisc.find("CoordinatesCentre")
            d["Mcircdisc_" + str(i)] = [eval(CoordinatesCentre.find("x").text), eval(CoordinatesCentre.find("y").text), eval(CoordinatesCentre.find("z").text)]
            d["Rcircdisc_" + str(i)] = eval(CircularDisc.find("Radius").text)
            Orientation = CircularDisc.find("Orientation")
            d["Phicircdisc_" + str(i)] = math.radians(eval(Orientation.find("Phi").text))
            d["Thetacircdisc_" + str(i)] = math.radians(eval(Orientation.find("Theta").text))
            d["Vcircdisc_" + str(i)] = eval(CircularDisc.find("Potential").text)

            i += 1

    i = 1
    for CircularDisc4Holes in root.iter("CircularDisc4Holes"):
        if CircularDisc4Holes.attrib["status"] == "enabled":
            CoordinatesCentre = CircularDisc4Holes.find("CoordinatesCentre")
            d["Mcircdisc4h_" + str(i)] = [eval(CoordinatesCentre.find("x").text), eval(CoordinatesCentre.find("y").text), eval(CoordinatesCentre.find("z").text)]
            d["Rcircdisc4h_" + str(i)] = eval(CircularDisc4Holes.find("Radius").text)
            Orientation = CircularDisc4Holes.find("Orientation")
            d["Ncircdisc4h_" + str(i)] = [eval(Orientation.find("x").text), eval(Orientation.find("y").text), eval(Orientation.find("z").text)]
            d["Vcircdisc4h_" + str(i)] = eval(CircularDisc4Holes.find("Potential").text)
            Holes = CircularDisc4Holes.find("Holes")
            d["Rcircdisc4h1_" + str(i)] = eval(Holes.find("Radius1").text)
            d["Rcircdisc4h2_" + str(i)] = eval(Holes.find("Radius2").text)
            d["Rcircdisc4h3_" + str(i)] = eval(Holes.find("Radius3").text)
            d["Rcircdisc4h4_" + str(i)] = eval(Holes.find("Radius4").text)
            CoordinatesCentre1 = Holes.find("CoordinatesCentre1")
            d["Mcircdisc4h1_" + str(i)] = [eval(CoordinatesCentre1.find("x").text),eval(CoordinatesCentre1.find("y").text),eval(CoordinatesCentre1.find("z").text)]
            CoordinatesCentre2 = Holes.find("CoordinatesCentre2")
            d["Mcircdisc4h2_" + str(i)] = [eval(CoordinatesCentre2.find("x").text),eval(CoordinatesCentre2.find("y").text),eval(CoordinatesCentre2.find("z").text)]
            CoordinatesCentre3 = Holes.find("CoordinatesCentre3")
            d["Mcircdisc4h3_" + str(i)] = [eval(CoordinatesCentre3.find("x").text),eval(CoordinatesCentre3.find("y").text),eval(CoordinatesCentre3.find("z").text)]
            CoordinatesCentre4 = Holes.find("CoordinatesCentre4")
            d["Mcircdisc4h4_" + str(i)] = [eval(CoordinatesCentre4.find("x").text),eval(CoordinatesCentre4.find("y").text),eval(CoordinatesCentre4.find("z").text)]

            i += 1

    i = 1
    for Surface in root.iter("Surface"):
        if Surface.attrib["status"] == "enabled":
            x, y = sy.symbols('x y')
            d["SurfEq_" + str(i)] = eval(Surface.find("Equation").text)
            d["Vsurface_" + str(i)] = eval(Surface.find("Potential").text)
            i += 1

    i = 1
    for Line in root.iter("Line"):
        if Line.attrib["status"] == "enabled":
            CoordinatesPoint1 = Line.find("CoordinatesPoint1")
            d["Pline_" + str(i)] = [eval(CoordinatesPoint1.find("x").text), eval(CoordinatesPoint1.find("y").text),
                                eval(CoordinatesPoint1.find("z").text)]
            CoordinatesPoint2 = Line.find("CoordinatesPoint2")
            d["Rline_" + str(i)] = [eval(CoordinatesPoint2.find("x").text), eval(CoordinatesPoint2.find("y").text),
                                eval(CoordinatesPoint2.find("z").text)]
            d["Qline_" + str(i)] = eval(Line.find("Charge").text)
            i += 1

    Particle = root.find("Particle")
    d["v"] = eval(Particle.find("Velocity").text)
    d["Theta1"] = eval(Particle.find("Theta1").text)
    d["Phi1"] = eval(Particle.find("Phi1").text)
    d["a"] = eval(Particle.find("Acceleration").text)
    d["Theta2"] = eval(Particle.find("Theta2").text)
    d["Phi2"] = eval(Particle.find("Phi2").text)

    Position = Particle.find("Position")
    d["Position"] = [eval(Position.find("x").text), eval(Position.find("y").text), eval(Position.find("z").text)]

    if Particle.find("Type").text == "custom":
        d["Mass"] = eval(Particle.find("Mass").text)
        d["Charge"] = eval(Particle.find("Charge").text)
    else:
        if Particle.find("Type").text in ParticleDictionary:
            d["Mass"] = ParticleDictionary[Particle.find("Type").text][0]
            d["Charge"] = ParticleDictionary[Particle.find("Type").text][1]
        else:
            print("Particle not found in dictionary")


    Setup = root.find("Setup")

    Trajectory = Setup.find("Trajectory")
    TrajectoryBoundaries = Trajectory.find("TrajectoryBoundaries")
    d["xmin"] = eval(TrajectoryBoundaries.find("xmin").text)
    d["xmax"] = eval(TrajectoryBoundaries.find("xmax").text)
    d["ymin"] = eval(TrajectoryBoundaries.find("ymin").text)
    d["ymax"] = eval(TrajectoryBoundaries.find("ymax").text)
    d["zmin"] = eval(TrajectoryBoundaries.find("zmin").text)
    d["zmax"] = eval(TrajectoryBoundaries.find("zmax").text)
    d["t"] = eval(Trajectory.find("TimeSteps").text)
    d["maxtime"] = eval(Trajectory.find("TimeLimit").text)

    Output = Setup.find("Output")
    d["TrajectoryPlot"] = Output.find("TrajectoryPlot").attrib["execute"]
    d["MagneticFieldPlot"] = Output.find("MagneticFieldPlot").attrib["execute"]
    d["NormalizeMagneticFieldPlot"] = Output.find("MagneticFieldPlot").attrib["normalize"]
    d["ElectricFieldPlot"] = Output.find("ElectricFieldPlot").attrib["execute"]
    d["NormalizeElectricFieldPlot"] = Output.find("ElectricFieldPlot").attrib["normalize"]
    d["WriteDataToFile"] = Output.find("WriteDataToFile").attrib["execute"]

    MagneticFieldPlot = Output.find("MagneticFieldPlot")
    MagneticFieldBoundaries =  MagneticFieldPlot.find("MagneticFieldBoundaries")
    d["xmin1"] = eval(MagneticFieldBoundaries.find("xmin").text)
    d["xmax1"] = eval(MagneticFieldBoundaries.find("xmax").text)
    d["ymin1"] = eval(MagneticFieldBoundaries.find("ymin").text)
    d["ymax1"] = eval(MagneticFieldBoundaries.find("ymax").text)
    d["zmin1"] = eval(MagneticFieldBoundaries.find("zmin").text)
    d["zmax1"] = eval(MagneticFieldBoundaries.find("zmax").text)

    ElectricFieldPlot = Output.find("ElectricFieldPlot")
    ElectricFieldBoundaries =  ElectricFieldPlot.find("ElectricFieldBoundaries")
    d["xmin2"] = eval(ElectricFieldBoundaries.find("xmin").text)
    d["xmax2"] = eval(ElectricFieldBoundaries.find("xmax").text)
    d["ymin2"] = eval(ElectricFieldBoundaries.find("ymin").text)
    d["ymax2"] = eval(ElectricFieldBoundaries.find("ymax").text)
    d["zmin2"] = eval(ElectricFieldBoundaries.find("zmin").text)
    d["zmax2"] = eval(ElectricFieldBoundaries.find("zmax").text)

    d["FileName"] = Output.find("WriteDataToFile").text

    return d


def SetupElements(d):
    i = 1
    Blist = []


    while "P_" + str(i) in d:
        Blist.append(StraightConductor(d["P_" + str(i)], d["Q_" + str(i)], d["I_" + str(i)]))
        i += 1


    i = 1
    j = 1

    while str(j) + "S_" + str(i) in d:
        list = ()
        while str(j) + "S_" + str(i) in d:
            list += (d[str(j) + "S_" + str(i)],)
            i += 1
        Blist.append(StraightConductorCollection(d["Icollection_" + str(j)], *list))
        i = 1
        j += 1


    i = 1
    while "SP_" + str(i) in d:
        Blist.append(Rectangularcoil(d["w_" + str(i)], d["l_" + str(i)], d["h_" + str(i)], d["N_" + str(i)], d["SP_" + str(i)], d["Phireccoil_" + str(i)],  d["Thetareccoil_" + str(i)], d["Psireccoil_" + str(i)], d["Ireccoil_" + str(i)], d["begin_" + str(i)]))
        i += 1

    i = 1

    while "M_" + str(i) in d:
        Blist.append(CircularConductor(d["M_" + str(i)], d["R_" + str(i)], d["Phi_" + str(i)], d["Theta_" + str(i)], d["Icircle_" + str(i)]))
        i += 1


    i = 1

    while "Mbent_" + str(i) in d:
        Blist.append(BentConductor(d["Mbent_" + str(i)], d["Rbent_" + str(i)], d["Phibent_" + str(i)], d["Thetabent_" + str(i)], d["Interval_" + str(i)], d["Ibent_" + str(i)]))
        i += 1


    i = 1

    while "Mcoil_" + str(i) in d:
        Blist.append(CircularCoil(d["Mcoil_" + str(i)], d["Rcoil_" + str(i)], d["Phicoil_" + str(i)], d["Thetacoil_" + str(i)], d["begincoil_" + str(i)], d["hcoil_" + str(i)], d["Ncoil_" + str(i)],  d["Icoil_" + str(i)]))
        i += 1


    return Blist


def SetupElements2(d):
    i = 1
    Elist = []


    while "Msphere_" + str(i) in d:
        Elist.append(Sphere(d["Msphere_" + str(i)], d["Rsphere_" + str(i)], d["Qsphere_" + str(i)]))
        i += 1


    i = 1
    while "Pline_" + str(i) in d:
        Elist.append(Line(d["Pline_" + str(i)], d["Rline_" + str(i)], d["Qline_" + str(i)]))
        i += 1


    return Elist

def SetupElementsWOS(d):
    ElistWOS = []

    i = 1
    while "Msphere2_" + str(i) in d:
        ElistWOS.append(Sphere2(d["Msphere2_" + str(i)], d["Rsphere2_" + str(i)], d["Vsphere2_" + str(i)]))
        i += 1

    i = 1
    while "Pcylinder_" + str(i) in d:
        ElistWOS.append(Cylinder(d["Pcylinder_" + str(i)], d["Qcylinder_" + str(i)], d["Rcylinder_" + str(i)], d["Vcylinder_" + str(i)]))
        i += 1

    i = 1
    while "Mcircdisc_" + str(i) in d:
        ElistWOS.append(CircularDisk(d["Mcircdisc_" + str(i)], d["Rcircdisc_" + str(i)], d["Phicircdisc_" + str(i)], d["Thetacircdisc_" + str(i)], d["Vcircdisc_" + str(i)]))
        i += 1

    i = 1
    while "Mcircdisc4h_" + str(i) in d:
        ElistWOS.append(CircularDisk4Holes(d["Mcircdisc4h_" + str(i)], d["Mcircdisc4h1_" + str(i)], d["Mcircdisc4h2_" + str(i)], d["Mcircdisc4h3_" + str(i)], d["Mcircdisc4h4_" + str(i)], d["Rcircdisc4h_" + str(i)], d["Rcircdisc4h1_" + str(i)], d["Rcircdisc4h2_" + str(i)], d["Rcircdisc4h3_" + str(i)], d["Rcircdisc4h4_" + str(i)], d["Ncircdisc4h_" + str(i)], d["Vcircdisc4h_" + str(i)]))
        i += 1

    i = 1
    while "SurfEq_" + str(i) in d:
        ElistWOS.append(Surface(d["SurfEq_" + str(i)], d["Vsurface_" + str(i)]))
        i += 1

    return ElistWOS

def ResultingField(list):
    F = Vector.zero

    for value in list:
       F += value

    return F

def ResultingString(list):
    F = ""
    if list != None:

        for value in list:
           F += value
    return F
