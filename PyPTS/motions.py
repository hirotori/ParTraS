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

_pypts._pypts_flib.delete_motion.argtypes = []
_pypts._pypts_flib.delete_motion.restype = None

def disable_motion():
    """
    delete current motion object from simulator.
    """
    _pypts._pypts_flib.delete_motion()

def droplet(dt:float, Re:float, rho_f:float, rho_p:float, RK_order:int,
            gravity:list[float]) -> None:
    """
    droplet motion. Particles are considered as uniform spheres of radius r and density `rho_p`, 
    and they move in a fluid of density `rho_f` and Reynolds number `Re`. 
    Particles are subjected to gravity and hydrodynamic drag force from the fluid. 
    The magnitude of the drag force is a function of the Reynolds number, 
    which is determined by the Reynolds number of the flow field `Re` and the velocity of the particle. 
    The motion of particles is discretized and integrated by the k-step low-storage Runge-Kutta method. 

    The user can adjust the number of steps of the Runge-Kutta method through `RK_order`. In the case of `RK_order=1`, 
    it is the same as the first-order Euler method. 
    Usually, `RK_order=4` is often used. This is the same as the Runge-Kutta method with 4 stages and 4th order.     
    When `RK_order`=**0**, a lower-order precision scheme (called here "Legacy scheme")is used instead of the Runge-Kutta method. 
    Legacy scheme is stable even for very small particles (e.g. radius <= 1e-6 [m]). 

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
    if RK_order < 0:
        raise ValueError("RK_order must be >= 0")
    if dt < 0:
        raise ValueError("dt must be > 0")
    if rho_f < 0:
        raise ValueError("rho_f must be > 0")
    if rho_p < 0:
        raise ValueError("rho_p must be > 0")
    
    _gravity = np.array(gravity, dtype=np.float64).reshape([3])

    _pypts._pypts_flib.assign_droplet_motion(
                                 ctypes.byref(ctypes.c_double(dt)),
                                 ctypes.byref(ctypes.c_double(Re)),
                                 ctypes.byref(ctypes.c_double(rho_f)),
                                 ctypes.byref(ctypes.c_double(rho_p)),
                                 ctypes.byref(ctypes.c_int(RK_order)),
                                 _gravity)


    