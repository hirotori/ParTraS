import pyvista as pv
import numpy as np
grid = pv.UnstructuredGrid("sax_flow.vtk")

xfor = np.loadtxt("conn.txt")
x = np.concatenate((grid.offset, grid.cell_connectivity))+1

for xi, xfi in zip(x, xfor):
    if xi != xfi:
        print(xi, xfi)

cell_id = 18913 - 1
ids = grid.cell_neighbors(cell_id, "faces")
print([idx+1 for idx in ids])

print(grid.get_cell(18912).point_ids)
for idx in ids:
    print(grid.get_cell(idx).point_ids)