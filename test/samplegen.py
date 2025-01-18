import numpy as np
from vtkmodules import all as vtk
from vtkmodules.util import numpy_support

class SimpleStructuredGrid():
    def __init__(self, imax:int, jmax:int, kmax:int) -> None:
        self.grid_ = vtk.vtkStructuredGrid()
        self.grid_.SetDimensions(imax,jmax,kmax)
        self.points_ = vtk.vtkPoints()

    def CreateStructuredGrid(self, point_array:np.array):
        self.points_.SetData(numpy_support.numpy_to_vtk(point_array))
        self.grid_.SetPoints(self.points_)

    def WriteOut(self, filename, binary = False, old_version=False, unstructured=False):

        if not unstructured:        
            writer_ = vtk.vtkStructuredGridWriter()
        else:
            writer_ = vtk.vtkUnstructuredGridWriter()
            mergeCells_ = vtk.vtkMergeCells()
            ugrid = vtk.vtkUnstructuredGrid()
            mergeCells_.SetUnstructuredGrid(ugrid)
            mergeCells_.SetTotalNumberOfDataSets(1)
            mergeCells_.SetTotalNumberOfCells(self.grid_.GetNumberOfCells())
            mergeCells_.SetTotalNumberOfPoints(self.grid_.GetNumberOfPoints())
            mergeCells_.MergeDataSet(self.grid_)
            ugrid = mergeCells_.GetUnstructuredGrid()
            if np.isnan(ugrid.GetNumberOfPoints()):
                ValueError("Merge Error.")


        writeout_core_((self.grid_ if not unstructured else ugrid), writer_, filename, binary, old_version)


def writeout_core_(grid_,writer_,filename,binary,oldversion):
    '''
    Note
    ----
    From vtkDataReader.h;
    "Currently VTK can write out two different versions of file format: files
    of VTK reader version 4.2 and previous; and VTK reader version 5.1 and
    later. This will likely change in the future. (Note: the major
    difference in the two formats is the way cell arrays are written out.)
    By default, Version 5.1 files are written out."
    '''
    writer_.SetInputData(grid_)
    writer_.SetFileName(filename)

    if binary :
        writer_.SetFileTypeToBinary()
    else:
        writer_.SetFileTypeToASCII()
        
    if oldversion:
        # https://vtk.org/doc/nightly/html/vtkDataWriter_8h_source.html
        # SetFileVersionはvtk9.1.0以降対応. 
        if vtk.VTK_MAJOR_VERSION >= 9 and vtk.VTK_MINOR_VERSION >= 1:
            writer_.SetFileVersion(writer_.VTK_LEGACY_READER_VERSION_4_2)
        else:
            raise ValueError("VTK version error :: Current version is {0}.{1}".format(vtk.VTK_MAJOR_VERSION, vtk.VTK_MINOR_VERSION))

    writer_.Write()


Lx = Ly = 10.0
Lz = 1.0
Nx = Ny = 2+1
Nz = 2+1

x = np.linspace(0.0, Lx, num=Nx)
y = np.linspace(0.0, Ly, num=Ny)
z = np.linspace(0.0, Lz, num=Nz)
rx, ry, rz = np.meshgrid(x, y, z, indexing="ij")
xyz = np.vstack((rx.T.reshape([-1]),ry.T.reshape([-1]),rz.T.reshape([-1]))).T

grid = SimpleStructuredGrid(Nx, Ny, Nz)
grid.CreateStructuredGrid(xyz)
grid.WriteOut("sample_small2x2x2.vtk", unstructured=True)
grid.WriteOut("sample_small2x2x2_old.vtk", unstructured=True, old_version=True)