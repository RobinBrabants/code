# contains funtion which plots the magnetic and electric fields which are able to be analytically evaluated


from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from sympy.vector import CoordSys3D
import sympy as sy
from numpy import sqrt
import math
from scipy.integrate import quad
import warnings


def Plotfield(Field, FieldType, d):
    from Lib.Functions import UpdateDictionary
    # funtion which plots the magnetic and electric fields which are able to be analytically evaluated

    print("Plotting the " + FieldType + " field...\r\n")

    warnings.filterwarnings("ignore")

    Phi2, t = sy.symbols('Phi2 t')

    pi = math.pi
    mu = 4 * pi * 10 ** (-7)

    factor = 1 / (((mu * 1) / (2 * pi * 40)))

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    L = CoordSys3D('L')

    if FieldType == "magnetic":
         x, y, z = np.meshgrid(np.arange(d["xmin1"], d["xmax1"], abs(d["xmax1"] - d["xmin1"]) / 10), np.arange(d["ymin1"], d["ymax1"], abs(d["ymax1"] - d["ymin1"]) / 10), np.arange(d["zmin1"], d["zmax1"], abs(d["zmax1"] - d["zmin1"]) / 10))
    else:
        x, y, z = np.meshgrid(np.arange(d["xmin2"], d["xmax2"], abs(d["xmax2"] - d["xmin2"]) / 10), np.arange(d["ymin2"], d["ymax2"], abs(d["ymax2"] - d["ymin2"]) / 10), np.arange(d["zmin2"], d["zmax2"], abs(d["zmax2"] - d["zmin2"]) / 10))

    Fieldcomponents = Field.components
    UpdateDictionary(Fieldcomponents)

    for basis in [L.i, L.j, L.k]:

        sum_integrals = np.zeros(x.shape)

        if str(Fieldcomponents[basis]).find("Integral(") != -1:
            while str(str(Fieldcomponents[basis])).find("Integral(") != -1:

                end1 = str(str(Fieldcomponents[basis])).find("Integral(")

                integrand = str(Fieldcomponents[basis])[str(Fieldcomponents[basis]).find("Integral(") + 9:str(Fieldcomponents[basis]).find(", ")]
                integrand = sy.sympify(integrand)
                integrand = eval(str(integrand))

                info1 = str(Fieldcomponents[basis])[str(Fieldcomponents[basis]).find(",") + 1:]
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

                Fieldcomponents[basis] = str(Fieldcomponents[basis])[:end1] + info2[begin2 + 3:]

            Fieldcomponents[basis] = eval(Fieldcomponents[basis] + "0")

            Fieldcomponents[basis] = np.add(Fieldcomponents[basis], sum_integrals)

        else:
            Fieldcomponents[basis] = eval(str(Fieldcomponents[basis]))

    u = Fieldcomponents[L.i]
    v = Fieldcomponents[L.j]
    w = Fieldcomponents[L.k]

    if d["NormalizeMagneticFieldPlot"] == "True":
        ax.quiver(x, y, z, u, v, w, length=3, arrow_length_ratio=0.4, pivot="middle", normalize=True)
        fig.suptitle(FieldType + " field plot, normalized", fontsize=14, fontweight='bold')
    else:
        ax.quiver(x, y, z, u, v, w, length=int(factor), arrow_length_ratio=0.4, pivot="middle")
        fig.suptitle(FieldType + " field plot", fontsize=14, fontweight='bold')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.ion()
    plt.show()

