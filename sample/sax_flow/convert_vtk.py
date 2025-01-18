import pyvista as pv
import numpy as np

interval = 100
for n in range(1000):
    digit = (n+1)*interval
    buff = np.loadtxt(f"particle_{digit}.pdata", skiprows=3)
    pos = buff[:,0:3]
    npart = len(pos)
    status = buff[:,11]

    # [1 0 1 1 1 2 ... 1 npart]
    cells = np.ones(2*npart, dtype=np.int32); cells[1::2] = np.arange(npart)
    cell_types = np.ones(npart, dtype=np.int32)
    grid = pv.UnstructuredGrid(cells, cell_types, pos)
    grid.cell_data.set_scalars(status, name="status")
    grid.save(f"particle_{digit}.vtk", binary=False)
