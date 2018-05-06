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
from sympy import sin, cos


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

    if FieldType == "Magnetic":
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
                        from sympy import sin, cos
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

    if d["Normalize" + FieldType + "FieldPlot"] == "True":
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


def PlotTrajectorys(particles):
    # funtion which plots all the trajectorys of the particles in a single figure

    print("Plotting all of the trajectorys in 1 figure...\r\n")

    warnings.filterwarnings("ignore")

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    for particle in particles:

        x = particle.Trajectory["x"]
        y = particle.Trajectory["y"]
        z = particle.Trajectory["z"]
        t = particle.Trajectory["t"]

        xs = np.array(x)
        ys = np.array(y)
        zs = np.array(z)
        ts = np.array(t)

        ax.plot(x, y, z, label="%s (%s)" % (particle.name, particle.Type))
        # ax.scatter maybe also good

        len_t = len(t) - 1
        ind_pos = []

        for i in range(0, 11):
            ind_pos.append(0 + i * (int(len_t / 10)))  # time labels

        xx = (xs[ind_pos])
        yy = (ys[ind_pos])
        zz = (zs[ind_pos])
        tt = (ts[ind_pos])

        for t, x, y, z in zip(tt, xx, yy, zz):
            label = '%s' % round(t,4)
            ax.text(x, y, z, label)
            ax.scatter(x, y, z, c="red")

        if "collision" in particle.Trajectory:
            electrode = particle.Trajectory["collision"]
            fig.suptitle('trajectory of ' + particle.name + ' (' + particle.Type + '), Particle collided with: ' + str(
                electrode.name) + '\r\n\r\n', fontsize=14, fontweight='bold')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.legend()

    fig.suptitle('trajectory of all particles', fontsize=14, fontweight='bold')

    ax.set_title('time in seconds', style='italic', fontsize=8, bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 7})

    plt.show(block=False)

