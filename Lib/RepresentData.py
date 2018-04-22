from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from sympy.vector import CoordSys3D
from Lib.Functions import *
import sympy as sy
from numpy import sqrt
import math
from scipy.integrate import quad
from sympy import sin, cos


def Plotfield(B, d):
    # funtion which plots the magnetic and electric fields which are able to be analytically evaluated

    print("Plotting the trajectory of the particle")

    print("Plotfield")

    Phi2, t = sy.symbols('Phi2 t')

    pi = math.pi
    mu = 4 * pi * 10 ** (-7)

    factor = 1 / (((mu * 1) / (2 * pi * 40)))

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    L = CoordSys3D('L')

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

