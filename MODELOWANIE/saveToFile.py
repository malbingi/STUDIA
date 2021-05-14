def dumpStateVTR(_LX, _LY, _F, _CELLS, string_name, _UX, _UY, _map, _Rho):
    file = open(string_name, "w+")
    #file.write(str(iter))
    
    file.write("# vtk DataFile Version 2.0\n")
    file.write("2D-ADE data file \n")
    file.write("ASCII\n")
    file.write("DATASET RECTILINEAR_GRID\n")
    file.write("DIMENSIONS " + str(_LX) + " " + str(_LY) + " 1 \n")
    file.write("X_COORDINATES " + str(_LX) + " int\n")
    for i in range(_LX):
        file.write(str(i) + " ")
    file.write("\n")
    file.write("Y_COORDINATES " + str(_LY) + " int\n")
    for i in range(_LY):
        file.write(str(i) + " ")
    file.write("\n")
    file.write("Z_COORDINATES 1 int\n")
    file.write("0\n")
    file.write("POINT_DATA " + str(_LX*_LY) + "\n")
    file.write("SCALARS density double\n")
    file.write("LOOKUP_TABLE default\n")
    for y in range(_LY):
        for x in range(_LX):
            """ phi = 0.0
            for f in range(_F):
                phi += _CELLS[x,y,f] """
            if _map[x,y] != 1:
                file.write(str(_Rho[x,y]) + "\n")
            else:
                file.write(str(0) +"\n")
            #file.write(str(_CELLS[x,y]) + "\n")
    file.write("\n")
    file.write("VECTORS velocity double\n")
    for y in range(_LY):
        for x in range(_LX):
            """ phi = 0.0
            for f in range(_F):
                phi += _CELLS[x,y,f] """
            if _map[x,y] != 1:
                file.write(str(_UX[x,y]) + " " + str(_UY[x,y]) + " 0" + "\n")
            else:
                file.write(str(0) + " " + str(0) + " " + str(0) +"\n")
    file.write("\n")
    
    """
    for x in range(_LX):
        file.write(str(x) + " " + str(_CELLS[x,0]) + "\n")
    """
    file.close()