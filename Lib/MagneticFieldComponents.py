import math
import sympy as sy
from math import sqrt, sin, cos, tan
from sympy import solve_poly_system
from sympy.vector import CoordSys3D, Vector


def StraightConductor(P, Q, I):
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


def StraightConductorCollection(I, *args):
    B = Vector.zero

    for i in range(0, len(args) - 1):
        B += StraightConductor(args[i], args[i+1], I)
    return B


def Rectangularcoil(w, l, h, N, SP, Phi, Theta, Psi, I, begin):
    from Lib.Functions import UpdateDictionary
    """width, length, heigth, number of windings, starting point of windings, angles of euler, begin is either width or length"""
    #from .Functions import *

    Phi = math.radians(Phi)
    Theta = math.radians(Theta)
    Psi = math.radians(Psi)

    L = CoordSys3D('L')

    e_1 = (cos(Phi)*cos(Psi) - sin(Phi)*sin(Psi)*cos(Theta))*L.i + (sin(Phi)*cos(Psi) + cos(Phi)*sin(Psi)*cos(Theta))*L.j + sin(Theta)*sin(Psi)*L.k
    e_2 = (-cos(Phi)*sin(Psi) - sin(Phi)*cos(Psi)*cos(Theta))*L.i + (-sin(Phi)*sin(Psi) + cos(Phi)*cos(Psi)*cos(Theta))*L.j + sin(Theta)*cos(Psi)*L.k
    e_3 = (sin(Phi)*sin(Theta))*L.i + (-cos(Phi)*sin(Theta))*L.j + cos(Theta)*L.k

    if begin == "l":
        a_1 = l*e_1
        a_2 = w*e_2
        a_3_1 = ((l/(2*w + 2*l))*(h/N))*e_3
        a_3_2 = ((w/(2*w + 2*l))*(h/N))*e_3
    elif begin == "w":
        a_1 = w*e_1
        a_2 = l*e_2
        a_3_1 = ((w/(2*w + 2*l))*(h/N))*e_3
        a_3_2 = ((l/(2*w + 2*l))*(h/N))*e_3

    a_1 = a_1.components
    UpdateDictionary(a_1)
    a_2 = a_2.components
    UpdateDictionary(a_2)
    a_3_1 = a_3_1.components
    UpdateDictionary(a_3_1)
    a_3_2 = a_3_2.components
    UpdateDictionary(a_3_2)


    new_point = SP
    windings = 0
    list = []
    list.extend(SP)

    while windings < N:
        new_point[0] += a_1[L.i] + a_3_1[L.i]
        new_point[1] += a_1[L.j] + a_3_1[L.j]
        new_point[2] += a_1[L.k] + a_3_1[L.k]
        list.extend(new_point)
        windings += 0.25

        if windings >= N:
            break

        new_point[0] += a_2[L.i] + a_3_2[L.i]
        new_point[1] += a_2[L.j] + a_3_2[L.j]
        new_point[2] += a_2[L.k] + a_3_2[L.k]
        list.extend(new_point)
        windings += 0.25

        if windings >= N:
            break

        new_point[0] += -a_1[L.i] + a_3_1[L.i]
        new_point[1] += -a_1[L.j] + a_3_1[L.j]
        new_point[2] += -a_1[L.k] + a_3_1[L.k]
        list.extend(new_point)
        windings += 0.25

        if windings >= N:
            break

        new_point[0] += -a_2[L.i] + a_3_2[L.i]
        new_point[1] += -a_2[L.j] + a_3_2[L.j]
        new_point[2] += -a_2[L.k] + a_3_2[L.k]
        list.extend(new_point)
        windings += 0.25

    B = Vector.zero

    for i in range(0, len(list) - 5, 3):
        B += StraightConductor([list[i], list[i+1], list[i+2]], [list[i+3], list[i+4], list[i+5]], I)

    return B


def BentConductor(M, R, Phi, Theta, interval, I):
    x1 = M[0]
    y1 = M[1]
    z1 = M[2]

    Phi = math.radians(Phi)
    Theta = math.radians(Theta)

    pi = math.pi
    mu = 4 * pi * 10 ** (-7)

    L = CoordSys3D('L')

    x, y, z, Phi2 = sy.symbols('x y z Phi2')

    a = R*(sy.sin(Phi2)*sy.sin(Phi) + sy.cos(Phi2)*sy.cos(Theta)*sy.cos(Phi))*L.i + R*(-sy.sin(Phi2)*sy.cos(Phi) + sy.cos(Phi2)*sy.cos(Theta)*sy.sin(Phi))*L.j + R*(-sy.cos(Phi2)*sy.sin(Theta))*L.k

    r = (x - x1 + R*sy.cos(Phi2)*sy.sin(Phi) - R*sy.sin(Phi2)*sy.cos(Theta)*sy.cos(Phi))*L.i + (y - y1 -R*sy.cos(Phi2)*sy.cos(Phi) - R*sy.sin(Phi2)*sy.cos(Theta)*sy.sin(Phi))*L.j + (z - z1 + R*sy.sin(Phi2)*sy.sin(Theta))*L.k

    b = ((mu*I)/(4*pi))*((a.cross(r))/(r.magnitude()**3))

    B = sy.integrate(-b, (Phi2, math.radians(interval[0]), math.radians(interval[1])))

    return B


def CircularConductor(M, R, Phi, Theta, I):
    x1 = M[0]
    y1 = M[1]
    z1 = M[2]

    Phi = math.radians(Phi)
    Theta = math.radians(Theta)

    pi = math.pi
    mu = 4 * pi * 10 ** (-7)

    L = CoordSys3D('L')

    x, y, z, Phi2 = sy.symbols('x y z Phi2')

    a = R*(sy.sin(Phi2)*sy.sin(Phi) + sy.cos(Phi2)*sy.cos(Theta)*sy.cos(Phi))*L.i + R*(-sy.sin(Phi2)*sy.cos(Phi) + sy.cos(Phi2)*sy.cos(Theta)*sy.sin(Phi))*L.j + R*(-sy.cos(Phi2)*sy.sin(Theta))*L.k

    r = (x - x1 + R*sy.cos(Phi2)*sy.sin(Phi) - R*sy.sin(Phi2)*sy.cos(Theta)*sy.cos(Phi))*L.i + (y - y1 -R*sy.cos(Phi2)*sy.cos(Phi) - R*sy.sin(Phi2)*sy.cos(Theta)*sy.sin(Phi))*L.j + (z - z1 + R*sy.sin(Phi2)*sy.sin(Theta))*L.k

    b = ((mu*I)/(4*pi))*((a.cross(r))/(r.magnitude()**3))

    B = sy.integrate(-b, (Phi2, 0, 2*pi))

    return B

def CircularCoil(M, R, Phi, Theta, begin, h, N, I):
    from Lib.Functions import UpdateDictionary, ConvertAnglesToVector, ConvertVectorToAngles
    Phi = math.radians(Phi)
    Theta = math.radians(Theta)

    pi = math.pi
    mu = 4 * pi * 10 ** (-7)

    h = 0.25*(h/N)

    L = CoordSys3D('L')

    n = ConvertAnglesToVector(Phi, Theta)

    M = M[0]*L.i + M[1]*L.j + M[2]*L.k
    M = M + h*n


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
        k = solve_poly_system([kx*nx + ky*ny, kx**2 + ky**2 - 1], kx, ky)
        kx = k[0][0]
        ky = k[0][1]
    else:
        kx = 0
        ky = 1

    tan_alpha = h/R

    dnx, dny, dnz = sy.symbols('dnx, dny, dnz')

    dn = solve_poly_system([nx*dnx + ny*dny + nz*dnz, dnx**2 + dny**2 + dnz**2 - tan_alpha**2, (kx/tan_alpha)*dnx + (ky/tan_alpha)*dny - math.cos(begin)], dnx, dny, dnz)

    if len(dn) >= 2:
        dnx = dn[1][0]
        dny = dn[1][1]
        dnz = dn[1][2]
    else:
        dnx = dn[0][0]
        dny = dn[0][1]
        dnz = dn[0][2]

    dn = dnx*L.i + dny*L.j + dnz*L.k

    r = n + dn

    r = r.normalize()

    Phi, Theta = ConvertVectorToAngles(r)

    alpha = math.atan(tan_alpha)

    R = R/math.cos(alpha)

    B = Vector.zero
    windings = 0

    while windings < N:
        B += BentConductor(M, R, Phi, Theta, [begin, begin + 180], I)

        M = M[0] * L.i + M[1] * L.j + M[2] * L.k
        M = M + 2*h * n
        Mcomponents = M.components
        UpdateDictionary(Mcomponents)
        M = [Mcomponents[L.i], Mcomponents[L.j], Mcomponents[L.k]]

        dn = -dn
        r = n + dn
        r = r.normalize()
        Phi, Theta = ConvertVectorToAngles(r)

        begin = begin + 180

        windings += 0.5

    return B