import numpy as np
import ctypes
from . import _pypts

_pypts._pypts_flib.update_field_vtk.argtypes = [ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.c_char_p,
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_bool),
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_bool)]
_pypts._pypts_flib.update_field_vtk.restype = None

_pypts._pypts_flib.no_update_field.argtypes = None
_pypts._pypts_flib.no_update_field.restype = None


class vtkUpdater:
    def __init__(self, dt_f:float, dt_p:float, basename:str, pad:int, field_only:bool, interval:int, ascii:bool) -> None:
        self.dt_f = dt_f
        self.dt_p = dt_p
        self.basename = basename
        self.pad = pad
        self.field_only = field_only
        self.interval = interval
        self.ascii = ascii

    @classmethod
    def no_update(cls):
        """
        disable updater without setting any parameters
        """
        cls.disable(cls)

    def disable(self):
        _pypts._pypts_flib.no_update_field()
        _pypts.logger.info("updater disabled")