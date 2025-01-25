import numpy as np
import ctypes
from . import _pypts

_pypts._pypts_flib.assign_droplet_motion.argtypes = [
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_int),
                                         np.ctypeslib.ndpointer(dtype=np.float64)]
_pypts._pypts_flib.assign_droplet_motion.restype = None



def droplet(dt:float, L_ref:float, U_ref:float, rho_f:float, mu_f:float, rho_p:float, RK_order:int,
            gravity:list=[0.0, 0.0, -9.806]) -> None:
    """
    droplet motion. 
    """
    _pypts._pypts_flib.assign_droplet_motion(
                                 ctypes.byref(ctypes.c_double(dt)),
                                 ctypes.byref(ctypes.c_double(L_ref)),
                                 ctypes.byref(ctypes.c_double(U_ref)),
                                 ctypes.byref(ctypes.c_double(rho_f)),
                                 ctypes.byref(ctypes.c_double(mu_f)),
                                 ctypes.byref(ctypes.c_double(rho_p)),
                                 ctypes.byref(ctypes.c_int(RK_order)),
                                 np.array(gravity, dtype=np.float64))


    