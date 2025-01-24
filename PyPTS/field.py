import numpy as np
import ctypes
from . import _pypts

_pypts._pypts_flib.init_field_vtk.argtypes = [ctypes.c_char_p,
                                              ctypes.POINTER(ctypes.c_bool)]
_pypts._pypts_flib.init_field_vtk.restype = None

_pypts._pypts_flib.init_field_afdet.argtypes = [ctypes.c_char_p]
_pypts._pypts_flib.init_field_afdet.restype = None

def init_field_vtk(filename:str, read_ascii:bool):

    b_filename = filename.encode()
    _pypts._pypts_flib.init_field_vtk(b_filename,
                                      ctypes.byref(ctypes.c_bool(read_ascii)))

    _pypts.logger.info("field data succesfully initialized from {0}".format(filename))

def init_field_afdet(filename:str):

    b_filename = filename.encode()
    _pypts._pypts_flib.init_field_afdet(b_filename)

    _pypts.logger.info("field data succesfully initialized from {0}".format(filename))
