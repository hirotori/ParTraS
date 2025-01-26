import numpy as np
import ctypes
from . import _pypts

_pypts._pypts_flib.assign_droplet_motion.argtypes = [
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_int),
                                         np.ctypeslib.ndpointer(dtype=np.float64)]
_pypts._pypts_flib.assign_droplet_motion.restype = None



def droplet(dt:float, Re:float, rho_f:float, rho_p:float, RK_order:int,
            gravity) -> None:
    """
    droplet motion. Particles are considered as uniform spheres of radius r and density `rho_p`, 
    and they move in a fluid of density `rho_f`. 
    Particles are subjected to gravity and hydrodynamic drag force from the fluid. The magnitude of the drag force is a function of the Reynolds number, 
    which is determined by the Reynolds number of the flow field `Re` and the velocity of the particle. 
    The motion of particles is discretized and integrated by the k-step low-storage Runge-Kutta method. 

    Note
    ------
    Note that all input parameters must be non-dimensionalized by the corresponding representative parameters. 
    PyPTS has no way to determine if the dimensions of these parameters are appropriate. 
    The user is solely responsible for the adequacy of the parameters specified. 
    
    Parameters
    -----------
    dt (float) : time stepping size [LU^-1]
    Re (float) : Reynolds number
    rho_f (float) : density of fluid [ρ]
    rho_p (float) : density of particle [ρ]
    RK_order (int) : order of the Runge-Kutta scheme
    gravity (list[float]) : gravity force exerted on a particle. 
    """
    _pypts._pypts_flib.assign_droplet_motion(
                                 ctypes.byref(ctypes.c_double(dt)),
                                 ctypes.byref(ctypes.c_double(Re)),
                                 ctypes.byref(ctypes.c_double(rho_f)),
                                 ctypes.byref(ctypes.c_double(rho_p)),
                                 ctypes.byref(ctypes.c_int(RK_order)),
                                 np.array(gravity, dtype=np.float64))


    