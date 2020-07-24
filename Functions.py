import numpy as np
import vtk




def scale(m, r_max, t_max):
     return int(np.round(((m / r_max) * t_max)))

def creatPolyData(name,points, cellArray):
    polyData = vtk.vtkPolyData()
    writer_Polydata = vtk.vtkXMLPolyDataWriter()
    fileName_Polydata = "/home/yaser/Desktop/ExtractBoundries/"+name+"polyData.vtp"
    polyData.SetPoints(points)
    polyData.SetVerts(cellArray)
    writer_Polydata = vtk.vtkXMLPolyDataWriter()
    writer_Polydata.SetFileName(fileName_Polydata)
    writer_Polydata.SetInputData(polyData)
    writer_Polydata.Write()


def createHistogram(fileName,length,lenU,lenV,lenW ,bin,points):

    filename = "/home/yaser/Desktop/ExtractBoundries/histo" + fileName + ".vti"
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(bin+1, bin+1,bin+1)
    imageData.AllocateScalars(vtk.VTK_FLOAT, 1)

    dims = imageData.GetDimensions()

    for z in range(dims[2]):
        for y in range(dims[1]):
            for x in range(dims[0]):
                imageData.SetScalarComponentFromFloat(x, y, z, 0, 0.0)





    dict = {}
    for j in range(length):
        i=points.GetPoint(j)
        i_0 = scale(i[0],lenU,bin)
        i_1 = scale(i[1],lenV,bin)
        i_2 = scale(i[2],lenW,bin)
        s = str(i_0) +" "+str(i_1) +" "+str(i_2)
        if dict.get(s):
            dict[s] = dict[s]+1

        else:
            imageData.SetScalarComponentFromDouble(i_0, i_1, i_2, 0, 1)
            # imageData.GetScalarComponentAsFloat(x, y, z, 0)
            dict[s] = 1
    dictValues = dict.values()
    sortedSet = sorted(set(dictValues))
    dictScale={}
    count =0
    for i in range(len(sortedSet)):
        dictScale[sortedSet[i]]=i
    for i in dict:
        #imageData.SetScalarComponentFromDouble(int(j[0]), int(j[1]), int(j[2]), 0, dict[i])
        j = i.split()
        # print("i: {},j: {},K: {}, dict: {}".format(int(j[0]), int(j[1]), int(j[2]), dictScale[dict[i]]+10))
        # imageData.SetScalarComponentFromFloat(int(j[0]), int(j[1]), int(j[2]), 0, float(dictScale[dict[i]]+10))
        imageData.SetScalarComponentFromFloat(int(j[0]), int(j[1]), int(j[2]), 0, float(dict[i]))

        count+=1







    # dims = imageData.GetDimensions()
    # for i in range(dims[0]):
    #     for j in range(dims[1]):
    #         for k in range(dims[2]):
    #             print("I: {}, j: {}, k: {}, v: {}".format(i, j, k, imageData.GetScalarComponentAsDouble(i, j, k, 0)))



    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(imageData)
    writer.Write()
    print(imageData.GetDimensions())





def projectTetrahedron(pointIds,color,volDic,unstructuredGrid,colors):

    # tetraPoints = vtk.vtkPoints()
    # for i in range(4):
    #     tetraPoints.InsertNextPoint(points[i][0],points[i][1],points[i][2])


    tetra = vtk.vtkTetra()
    for i in range (4):
        tetra.GetPointIds().SetId(i, pointIds[i])

    cellId = unstructuredGrid.InsertNextCell(tetra.GetCellType(),tetra.GetPointIds())
    r = volDic[color]
    if color != 0:
        if r >1000:
            colors.InsertTuple1(cellId, 999)
        else:
            colors.InsertTuple1(cellId, volDic[color])


    else:
        colors.InsertTuple1(cellId, 0)





def computeVolumeOfTetrahedron(cellPointsValue):
    vertexA,vertexB,vertexC,vertexD = cellPointsValue[0],cellPointsValue[1],cellPointsValue[2],cellPointsValue[3]
    AB = [vertexB[0]-vertexA[0],vertexB[1]-vertexA[1],vertexB[2]-vertexA[2]]
    AC =  [vertexC[0]-vertexA[0],vertexC[1]-vertexA[1],vertexC[2]-vertexA[2]]
    AD = [vertexD[0]-vertexA[0],vertexD[1]-vertexA[1],vertexD[2]-vertexA[2]]
    vp= np.dot(AD,np.cross(AB,AC))
    vol =np.abs(vp)/6
    return vol







def process(cellPointsValue,pointIds,volDic,ren,unstructuredGrid, colors):
    vol = computeVolumeOfTetrahedron(cellPointsValue)
    if (vol!=0):
        projectTetrahedron(pointIds, vol,volDic,unstructuredGrid,colors)

    else:
        projectTetrahedron( pointIds, 0,volDic,unstructuredGrid,colors)





def pointInTriangle(v0,v1,v2,p):
    u = v1 - v0
    v = v2 - v0
    n = np.cross(u,v)
    w = p- v0
    npCross =np.cross(u,w)
    npDot =np.dot(n, n)
    gamma = np.dot(npCross,n)/npDot
    npCross= np.cross(w,v)
    npDot =np.dot(n,n)
    beta = np.dot(npCross,n)/npDot
    alpha = 1 - gamma - beta
    return ((0 <= alpha) and (alpha <= 1) and (0 <= beta) and (beta <= 1) and (0 <= gamma) and (gamma <= 1))




def sameside(v1,v2,v3,v4,p):
    normal = np.cross(v2 - v1, v3 - v1)
    return (np.dot(normal, v4 - v1) * np.dot(normal, p - v1) > 0)

def tetraCoord(A,B,C,D):
    v1 = B-A ; v2 = C-A ; v3 = D-A
    # mat defines an affine transform from the tetrahedron to the orthogonal system
    mat = np.array((v1,v2,v3)).T
    # The inverse matrix does the opposite (from orthogonal to tetrahedron)
    M1=mat
    try:
        M1 = np.linalg.inv(mat)
    except :
        print("error Occured")

    return(M1)

def pointInsideT(v1,v2,v3,v4,p):
    # # Find the transform matrix from orthogonal to tetrahedron system
    # M1 = tetraCoord(v1, v2, v3, v4)
    # # apply the transform to P
    # newp = M1.dot(p - v1)
    # # perform test
    # return (np.all(newp >= 0) and np.all(newp <= 1) and np.sum(newp) <= 1)
    return sameside(v1, v2, v3, v4, p) and sameside(v2, v3, v4, v1, p) and sameside(v3, v4, v1, v2, p) and sameside(v4, v1, v2, v3, p)

def incrmentCell(j,dic):

    s = str(j[0]) + " " + str(j[1]) + " " + str(j[2])
    if (s) in dic:
        co = dic[s] + 1
        dic[s] = co

    else:
        dic[s] = 1

def Bresenham3D(x1, y1, z1, x2, y2, z2):
    ListOfPoints = []
    ListOfPoints.append((x1, y1, z1))
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    dz = abs(z2 - z1)
    if (x2 > x1):
        xs = 1
    else:
        xs = -1
    if (y2 > y1):
        ys = 1
    else:
        ys = -1
    if (z2 > z1):
        zs = 1
    else:
        zs = -1

    # Driving axis is X-axis"
    if (dx >= dy and dx >= dz):
        p1 = 2 * dy - dx
        p2 = 2 * dz - dx
        while (x1 != x2):
            x1 += xs
            if (p1 >= 0):
                y1 += ys
                p1 -= 2 * dx
            if (p2 >= 0):
                z1 += zs
                p2 -= 2 * dx
            p1 += 2 * dy
            p2 += 2 * dz
            ListOfPoints.append((x1, y1, z1))

            # Driving axis is Y-axis"
    elif (dy >= dx and dy >= dz):
        p1 = 2 * dx - dy
        p2 = 2 * dz - dy
        while (y1 != y2):
            y1 += ys
            if (p1 >= 0):
                x1 += xs
                p1 -= 2 * dy
            if (p2 >= 0):
                z1 += zs
                p2 -= 2 * dy
            p1 += 2 * dx
            p2 += 2 * dz
            ListOfPoints.append((x1, y1, z1))

            # Driving axis is Z-axis"
    else:
        p1 = 2 * dy - dz
        p2 = 2 * dx - dz
        while (z1 != z2):
            z1 += zs
            if (p1 >= 0):
                y1 += ys
                p1 -= 2 * dz
            if (p2 >= 0):
                x1 += xs
                p2 -= 2 * dz
            p1 += 2 * dy
            p2 += 2 * dx
            ListOfPoints.append((x1, y1, z1))
    return ListOfPoints


def createHistogramFromUnstructuredGrid(f,data,bin,lenU,lenV,lenW,u,v,w,dictU,dictV,dictW):
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(bin + 1, bin + 1, bin + 1)
    imageData.AllocateScalars(vtk.VTK_DOUBLE, 1)
    filename = "/home/yaser/Desktop/ExtractBoundries/ExampleConinousResol"+f+".vti"

    dims = imageData.GetDimensions()

    for z in range(dims[2]):
        for y in range(dims[1]):
            for x in range(dims[0]):
                imageData.SetScalarComponentFromFloat(x, y, z, 0, 0.0)



    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(filename)



    dict = {}
    dict2={}
    h = 0
    print("in others...")
    for i in range(int(data.GetNumberOfCells())):
        h += 1
        cell = data.GetCell(i)
        cellPointsValue = []
        f2 = []
        for i_d in range(4):
            id = cell.GetPointId(i_d)
            j, k, l = dictU[u[id]], dictV[v[id]], dictW[w[id]]
            j = scale(j, lenU, bin)
            k= scale(k, lenV, bin)
            l= scale(l, lenW, bin)
            cellPointsValue.append((j, k, l))
        se = set(cellPointsValue)
        se1 = set()
        if len(se) == 1:
            p = se.pop()
            se1.add(p)
        elif len(se) == 2:


            p1 = se.pop()
            p2 = se.pop()
            se1.add(p1)
            se1.add(p2)
            l = Bresenham3D(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2])
            # for j in l:
            #     se1.add(j)
        elif len(se) == 3:



            p1 = se.pop()
            p2 = se.pop()
            p3 = se.pop()
            se1.add(p1)
            se1.add(p2)
            se1.add(p3)
            # l = Bresenham3D(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2])
            # for j in l:
            #     se1.add(j)
            # l = Bresenham3D(p1[0], p1[1], p1[2], p3[0], p3[1], p3[2])
            # for j in l:
            #     se1.add(j)
            # l = Bresenham3D(p3[0], p3[1], p3[2], p2[0], p2[1], p2[2])
            # for j in l:
            #     se1.add(j)
            l = [p1, p2, p3]
            lx = [x[0] for x in l]
            ly = [x[1] for x in l]
            lz = [x[2] for x in l]
            # print("minX: {}, maxX:{}, minY: {}, maxY:{}, minZ:{}, maxZ:{}".format(min(lx),max(lx),min(ly),max(ly),min(lz),max(lz)))
            for i in range(min(lx), max(lx) + 1):
                for j in range(min(ly), max(ly) + 1):
                    for k in range(min(lz), max(lz) + 1):
                        if pointInTriangle(np.array(p1), np.array(p2), np.array(p3), np.array([i, j, k])):
                            se1.add((i, j, k))
        else:

            p1 = se.pop()
            p2 = se.pop()
            p3 = se.pop()
            p4 = se.pop()
            # l = Bresenham3D(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2])
            # for j in l:
            #     se1.add(j)
            # l = Bresenham3D(p1[0], p1[1], p1[2], p3[0], p3[1], p3[2])
            # for j in l:
            #     se1.add(j)
            # l = Bresenham3D(p3[0], p3[1], p3[2], p2[0], p2[1], p2[2])
            # for j in l:
            #     se1.add(j)
            #
            # l = Bresenham3D(p1[0], p1[1], p1[2], p4[0], p4[1], p4[2])
            # for j in l:
            #     se1.add(j)
            # l = Bresenham3D(p4[0], p4[1], p4[2], p3[0], p3[1], p3[2])
            # for j in l:
            #     se1.add(j)
            # l = Bresenham3D(p4[0], p4[1], p4[2], p2[0], p2[1], p2[2])
            # for j in l:
            #     se1.add(j)
            l = [p1, p2, p3, p4]
            # print(l)
            lx = [x[0] for x in l]
            ly = [x[1] for x in l]
            lz = [x[2] for x in l]
            # print("minX: {}, maxX:{}, minY: {}, maxY:{}, minZ:{}, maxZ:{}".format(min(lx),max(lx),min(ly),max(ly),min(lz),max(lz)))
            for i in range(min(lx), max(lx) + 1):
                for j in range(min(ly), max(ly) + 1):
                    for k in range(min(lz), max(lz) + 1):
                        if pointInsideT(np.array(p1), np.array(p2), np.array(p3), np.array(p4), np.array([i, j, k])):
                            se1.add((i, j, k))


        while len(se1) > 0:
            p = se1.pop()
            incrmentCell(p, dict)
            if dict2.get(p):
                for q in f2:
                    dict2[p].append(q)

            else:
                dict2[p] = []
                for q in f2:
                    dict2[p].append(q)


        # print("h is: {}".format(h))
    values = dict.values()
    sortedSet = sorted(set(values))
    print("sortedSet {}".format(sortedSet))
    dict3 ={}

    for i in dict:
        # indx = int(sortedSet.index(dict[i])) + 1
        indx = int(dict[i]) + 1

        j = i.split(" ")
        a,b,c = int(j[0]), int(j[1]), int(j[2])
        imageData.SetScalarComponentFromDouble( a,b,c,0, indx)

        if dict3.get(indx):
            dict3[indx].extend(dict2[(a,b,c)])

        else:
            dict3[indx ] =[]
            dict3[indx].extend(dict2[(a,b,c)])





    writer.SetInputData(imageData)
    writer.Write()
    return imageData,dict3
