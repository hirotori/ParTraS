import numpy as np
import ctypes
from . import _pypts

_pypts._pypts_flib.assign_droplet_motion.argtypes = [ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_int)]
_pypts._pypts_flib.assign_droplet_motion.restype = None



def droplet(rho_f:float, mu_f:float, rho_p:float, dt:float, RK_order:int) -> None:

    _pypts._pypts_flib.assign_droplet_motion(ctypes.byref(ctypes.c_double(rho_f)),
                                 ctypes.byref(ctypes.c_double(mu_f)),
                                 ctypes.byref(ctypes.c_double(rho_p)),
                                 ctypes.byref(ctypes.c_double(dt)),
                                 ctypes.byref(ctypes.c_int(RK_order)))


    