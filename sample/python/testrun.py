from PyPTS import particle
from PyPTS import field
from PyPTS import motions
from PyPTS import update
from PyPTS import simulation
import numpy as np

#particle data
pdata = particle.ParticleData(10)

rng = np.random.default_rng(seed=1000)
r0 = [0.09, 0.2, 0.4]
pos = r0 + rng.random((10,3))*0.02 - 0.01
pdata.pos = pos
pdata.vel = np.zeros((10,3))
pdata.force = np.zeros((10,3))
pdata.radius = np.full(10, fill_value=1e-5)
pdata.ref_cell = np.ones(10)
pdata.state = np.ones(10)
pdata.initialize_particle_data()

#field data
field.init_field_vtk("../sax_flow/sax_flow.vtk", True)


# motion
n_rk = 4
rho_p = 1000.0
mu_p  = 0.001
rho_f = 1.0
mu_f  = 1e-5
dt = 0.0001
motions.droplet(rho_f, mu_f, rho_p, dt, n_rk)

# updater is disabled (i.e. steady flow)
update.vtkUpdater.no_update()

# simulation
simulation.initialize(100, True, "data", 
                      100, False, "data")

simulation.run(1, 1000)
