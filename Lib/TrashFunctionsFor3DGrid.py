



def WOS_Distance(P, functions, potentials):

    max_i = len(functions)

    dist = []
    #bepalen van de afstand vanuit P tot aan alle oppervlakken
    for i in range(0, max_i):
        if functions[i].find("sphere") != -1:
            M = eval(functions[i][functions[i].find("["):functions[i].find("]") + 1])
            R = eval(functions[i][functions[i].find(";") + 2:])
            dist.append(dist_sphere(P, M, R))


        elif functions[i].find("cylinder") != -1:
            P1 = eval(functions[i][functions[i].find("["):functions[i].find("]") + 1])
            Q = eval(functions[i][functions[i].find(";") + 2:functions[i].find(", radius")])
            R = eval(functions[i][functions[i].find("radius") + 9:])

            dist.append(dist_cylinder(P, P1, Q, R))

        elif functions[i].find("circulardisk") != -1:
            M = eval(functions[i][functions[i].find("["):functions[i].find("]") + 1])
            R = eval(functions[i][functions[i].find(";") + 2:functions[i].find(", orientation")])
            Phi = eval(functions[i][functions[i].find("phi") + 4:functions[i].find(", theta")])
            Theta = eval(functions[i][functions[i].find("theta") + 6:])

            dist.append(dist_circdisc(P, M, R, Phi, Theta))

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

            dist.append(dist_circdisc_4holes(P, M, M1, M2, M3, M4, R, R1, R2, R3, R4, N))

        else:
            function = eval(functions[i])
            dist.append(dist_surface(P, function))

    index, radius = min(enumerate(dist), key=itemgetter(1))

    return (potentials[index], radius)


def WOS_DistanceGrid(EWOS, dim_x, dim_y, dim_z):

    functions, potentials = Get_Data_EWOS(EWOS)

    Dgrid = np.empty(shape=(len(dim_z), len(dim_x), len(dim_y)), dtype=tuple)

    # multi-processing of the different z arrays in Dgrid (determines automatically how many parallel processes it can run depending on the amount of CPU's in your computer)
    pool = mp.Pool()
    results = [pool.apply_async(WOS_Process, args=(z, functions, potentials, dim_x, dim_y, dim_z)) for z in dim_z]
    pool.close
    pool.join

    list = [p.get() for p in results]

    # list has to be sorted cause processes don't finish in the correct order
    list.sort()
    list = [r[1] for r in list]

    i = 0
    for array in list:
        Dgrid[i,:,:] = array
        i+=1


    print("DistanceGrid: completed")

    np.set_printoptions(threshold=np.nan, linewidth=500)

    print(Dgrid.tolist())

    f = open("Dgrid", "w+")

    f.write("Dgrid\r\n")
    f.write(str(Dgrid.tolist()))


    return Dgrid.tolist()


def WOS_Process(z, functions, potentials, dim_x, dim_y, dim_z):

    Grid = np.empty(shape=(len(dim_x), len(dim_y)), dtype=tuple)

    print("DistanceGrid:", 100 * (dim_z.index(z) + 1) / len(dim_z), "%")

    for k in range(0, len(dim_y)):
        for l in range(0, len(dim_x)):
            Grid[k, l] = WOS_Distance([dim_x[l], dim_y[k], z], functions, potentials)

    return (z, Grid)


def WalkOnSpheres_potential_slice_DGrid(Dgrid, dim_x, dim_y, dim_z, space, x, y, z):

    Vgrid = np.zeros(shape=(len(z), len(y)))

    factor = len(y)
    number_iterations = len(y) * len(z)

    #multi-processing of the different points in Vgrid (determines automatically how many parallel processes it can run depending on the amount of CPU's in your computer)
    pool = mp.Pool()
    results = [pool.apply_async(MultiProcess_WOS_DGrid, args=([x, j, i], z[::-1].index(i), y.index(j), Dgrid, dim_x, dim_y, dim_z, space, factor, number_iterations)) for i in z[::-1] for j in y]
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


def MultiProcess_WOS_DGrid(P, k, l, Dgrid, dim_x, dim_y, dim_z, space, factor, number_iteration):
    potential = potential_EWOS_DGrid(Dgrid, dim_x, dim_y, dim_z, space, P)

    #to keep track of the process
    print("WOS process:",l + factor * k+1,"out of", number_iteration,"(",100*(l + factor * k+1)/number_iteration,"percent)")

    return (l + factor * k, potential)


def potential_EWOS_DGrid(Dgrid, dim_x, dim_y, dim_z, space, P):

    V_inf = 0
    maxsteps = 400  # indien er over dit aantal iteraties wordt gegaan, gooien we de iteratie weg.
    iterations = 200  # niet te verwarren met maxsteps

    V = 0
    k = 0
    values = 0

    while k < iterations:
        j = 0
        P_var = P[:]
        while j < maxsteps:
            P_var[0] = myround(P_var[0], space, prec=3)
            P_var[1] = myround(P_var[1], space, prec=3)
            P_var[2] = myround(P_var[2], space, prec=3)

            m = int((P_var[0]-dim_x[0])/space)
            n = int((P_var[1] - dim_y[0]) / space)
            o = int((P_var[2] - dim_z[0]) / space)

            #als de straal te groot is wordt de ptentiaal gelijk gesteld aan V_inf
            if o>len(Dgrid)-1 or o<0 or n>len(Dgrid[0])-1 or n<0 or m>len(Dgrid[0][0])-1 or m<0:
                V += V_inf

                values += 1

                break

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


radius = Dgrid[o][m][n][1]

# als het de eerste keer kleiner is dan space hebben we in feite de particlemove functie al stopgezet
if radius <= space:
    potential = Dgrid[o][m][n][0]

    V += potential

    values += 1

    break

else:
    # random vector aanmaken en P_var updaten
    x = random() - 0.5
    y = random() - 0.5
    z = random() - 0.5


def myround(x, base, prec=3):
  return round(base * round(float(x)/base),prec)