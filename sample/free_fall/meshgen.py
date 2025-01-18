import pyvista as pv
import numpy as np

Nx = Ny = 1
Nz = 100
x = np.linspace(0.0, 1.0, Nx+1)
y = np.linspace(0.0, 1.0, Ny+1)
z = np.linspace(0.0, 100.0, Nz+1)
v = np.zeros((Nx*Ny*Nz,3))

rx, ry, rz = np.meshgrid(x, y, z, indexing="ij")

grid = pv.StructuredGrid(rx, ry, rz)
ugrid = pv.UnstructuredGrid(grid)
ugrid.cell_data.set_vectors(v, name="velocity")
ugrid.save("quiescent.vtk", binary=False)
