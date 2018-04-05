from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from sympy.vector import CoordSysCartesian
from Lib.Functions import *
import sympy as sy
from numpy import sqrt
import math
from scipy.integrate import quad
from sympy import sin, cos


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

    for i in range(0,11):
        ind_pos.append(0 + i*(int(len_t/10)))

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
    ax.set_title('time in seconds', style='italic', fontsize=8, bbox={'facecolor':'red', 'alpha':0.2, 'pad':7})

    plt.show(block="False")


def Plotfield(B, d):

    print("Plotfield")

    Phi2, t = sy.symbols('Phi2 t')

    pi = math.pi
    mu = 4 * pi * 10 ** (-7)

    factor = 1 / (((mu * 1) / (2 * pi * 40)))

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    L = CoordSysCartesian('L')

    x, y, z = np.meshgrid(np.arange(d["xmin1"], d["xmax1"], abs(d["xmax1"] - d["xmin1"]) / 10),
                          np.arange(d["ymin1"], d["ymax1"], abs(d["ymax1"] - d["ymin1"]) / 10),
                          np.arange(d["zmin1"], d["zmax1"], abs(d["zmax1"] - d["zmin1"]) / 10))

    Bcomponents = B.components
    UpdateDictionary(Bcomponents)

    for basis in [L.i, L.j, L.k]:

        sum_integrals = np.zeros(x.shape)

        if str(Bcomponents[basis]).find("Integral(") != -1:
            while str(str(Bcomponents[basis])).find("Integral(") != -1:

                end1 = str(str(Bcomponents[basis])).find("Integral(")

                integrand = str(Bcomponents[basis])[str(Bcomponents[basis]).find("Integral(") + 9:str(Bcomponents[basis]).find(", ")]
                integrand = sy.sympify(integrand)
                integrand = eval(str(integrand))

                info1 = str(Bcomponents[basis])[str(Bcomponents[basis]).find(",") + 1:]
                dx = info1[info1.find("(") + 1:info1.find(",")]
                info2 = info1[info1.find(",") + 2:]
                a = info2[:info2.find(",")]
                b = info2[info2.find(",") + 2:info2.find("))")]

                begin2 = info2.find("))")

                print("caculating...")

                if dx == "Phi2":
                    for index, w in np.ndenumerate(integrand):
                        sum_integrals[index] += quad(lambda Phi2: eval(str(w)), eval(a), eval(b))[0]
                if dx == "t":
                    for index, w in np.ndenumerate(integrand):
                        sum_integrals[index] += quad(lambda t: eval(str(w)), eval(a), eval(b))[0]

                Bcomponents[basis] = str(Bcomponents[basis])[:end1] + info2[begin2 + 3:]

            print("done")

            Bcomponents[basis] = eval(Bcomponents[basis] + "0")

            Bcomponents[basis] = np.add(Bcomponents[basis], sum_integrals)

        else:
            Bcomponents[basis] = eval(str(Bcomponents[basis]))


    print("almost")
    u = Bcomponents[L.i]
    v = Bcomponents[L.j]
    w = Bcomponents[L.k]

    if d["NormalizeMagneticFieldPlot"] == "True":
        ax.quiver(x, y, z, u, v, w, length=3, arrow_length_ratio=0.4, pivot="middle", normalize=True)
    else:
        ax.quiver(x, y, z, u, v, w, length=int(factor), arrow_length_ratio=0.4, pivot="middle")

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    fig.suptitle('field', fontsize=14, fontweight='bold')

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