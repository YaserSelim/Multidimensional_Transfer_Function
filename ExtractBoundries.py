
import vtk
import numpy as np
import numpy.linalg as la
import Functions as fun
import sys
import random

name='20bbbb'
bin =50
points = vtk.vtkPoints()
unstructuredGrid = vtk.vtkUnstructuredGrid()


#VTK_POLYDATA

cellArray = vtk.vtkCellArray()

#End_VTK_POLYDATA



#VTK_IMAGEDATA
reader = vtk.vtkStructuredPointsReader()
reader.SetFileName("/home/yaser/Desktop/ExtractBoundries/dataset2.vtk")
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()
reader.Update()
data = reader.GetOutput()
# data.GetPointData().GetArray('data').SetName('vectors')
pointData = np.array(data.GetPointData().GetArray('vectors'))
dimX= data.GetDimensions()[0]






data2 = reader.GetOutput()
# data2 = vtk.vtkStructuredPoints()
# data2.SetDimensions(data.GetDimensions())
data2.AllocateScalars(vtk.VTK_DOUBLE, 1)


print(data.GetPointData().GetArray('vectors').GetTuple(6))
#print(data.GetScalarComponentAsDouble(20,3,17,0))
pointData = np.array(data.GetPointData().GetArray('vectors'))
gradient = vtk.vtkGradientFilter()
gradient.SetInputData(data)
gradient.SetInputArrayToProcess(0, 0, 0, 0, "vectors")
gradient.SetResultArrayName("grad")
gradient.Update()
g = gradient.GetOutput()
gradientValues =np.array(g.GetPointData().GetArray('grad'))
print(gradientValues[4])
print("pointdata: {}".format(data.GetPoint(61)))
unstructuredGrid.SetPoints(points)
colors = vtk.vtkFloatArray()
u,v,w = [],[],[]

for i in range(int(gradientValues.shape[0])):
    x = np.array(gradientValues[i].reshape((3, 3)))
    y = la.norm(x, 2)
    u0,y0,w0 =pointData[i][0],pointData[i][1],y
    u.append(u0)
    v.append(y0)
    w.append(w0)


// this way to get the pointData faster
uSortedSet = sorted(set(u))
vSortedSet = sorted(set(v))
wSortedSet = sorted(set(w))
print("len U: {}".format((wSortedSet)))
leng = len(u)
dictU={}
dictV={}
dictW={}
leng = len(u)
lenU,lenV,lenW = len(uSortedSet),len(vSortedSet),len(wSortedSet)
m = max(lenU,lenV,lenW)
for i in range(m):
    if i <lenU:
        dictU[uSortedSet[i]] = i
    if i<lenV:
        dictV[vSortedSet[i]]=i
    if i < lenW:
        dictW[wSortedSet[i]]=i



for i in range(leng):
    j, k, l = dictU[u[i]], dictV[v[i]], dictW[w[i]]
    id = points.InsertNextPoint(j, k, l)
    r = cellArray.InsertNextCell(1)
    r = cellArray.InsertCellPoint(id)


fun.creatPolyData(name,points,cellArray)
values = sorted(dictU.values())

fun.createHistogram(name,len(u),lenU,lenV,lenW,bin,points)

# sys.exit()

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(300, 150)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
ren.SetBackground(1, 1, 1)

volume=[]
print(len(u))
for i in range (int(data.GetNumberOfCells() )):
    cell = data.GetCell(i)
    cellPointsValue=[]
    pointIds =[]
    for i_d in range(4):
        id = cell.GetPointId(i_d)
        pointIds.append(id)
        j, k, l = dictU[u[id]],dictV[v[id]], dictW[w[id]]
        cellPointsValue.append([j,k,l])
    volume.append(fun.computeVolumeOfTetrahedron(cellPointsValue))


volume = sorted(set(volume),reverse=True)
volDic ={}
for i in range (len(volume)):
    volDic[volume[i]] = i


for i in range (int(data.GetNumberOfCells() )):
    cell = data.GetCell(i)
    cellPointsValue=[]
    pointIds =[]
    for i_d in range(4):
        id = cell.GetPointId(i_d)
        pointIds.append(id)
        j, k, l =dictU[u[id]],dictV[v[id]], dictW[w[id]]
        cellPointsValue.append([j,k,l])
    fun.process(cellPointsValue,pointIds,volDic,ren,unstructuredGrid,colors)
unstructuredGrid.GetCellData().SetScalars(colors)



imagedata2 ,dict2= fun.createHistogramFromUnstructuredGrid(name,data,bin, lenU,lenV,lenW,u,v,w,dictU,dictV,dictW)

print("End of prog")



# lut.SetTableValue(1, 0, 0, 1, .5) # Red
# lut.SetTableValue(2, 0, 1, 1, .5) # Green
# lut.SetTableValue(3, 0, 0, 1, .5) # Blue
# lut.SetTableValue(4, 1, 0, 0, .5) # White









aTetraMapper = vtk.vtkDataSetMapper()
aTetraMapper.SetScalarModeToUseCellData()
aTetraMapper.UseLookupTableScalarRangeOn()
aTetraMapper.SetInputData(unstructuredGrid)
writer =vtk.vtkXMLUnstructuredGridWriter()
fileName = "/home/yaser/Desktop/ExtractBoundries/unstr"+name+".vtu"
writer.SetFileName(fileName)
writer.SetInputData(unstructuredGrid)
writer.Write()

aTetraActor = vtk.vtkActor()
aTetraActor.SetMapper(aTetraMapper)
ren.AddActor(aTetraActor)
ren.ResetCamera()
iren.Initialize()
renWin.Render()
ren.GetActiveCamera().Zoom(1.5)
renWin.Render()
iren.Start()

