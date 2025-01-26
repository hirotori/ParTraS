import numpy as np
import ctypes
from . import _pypts

_pypts._pypts_flib.init_field_vtk.argtypes = [ctypes.c_char_p,
                                              ctypes.POINTER(ctypes.c_bool)]
_pypts._pypts_flib.init_field_vtk.restype = None

_pypts._pypts_flib.init_field_afdet.argtypes = [ctypes.c_char_p]
_pypts._pypts_flib.init_field_afdet.restype = None

def init_field_vtk(filename:str, read_ascii:bool):
    """
    Initializes the flow field. 

    The flow field is constructed from grid data stored in Legacy VTK unstructured grid format (.vtk). 
    If the grid data also contains field data such as velocity and temperature fields, they are restored as well. 

    Parameters
    -----------
    filename (str) : filename in which the grid data is stored.
    read_ascii (bool) : read the file in ASCII format.

    Note
    ------
    The name of the field data stored in the flow field must be `velocity` for velocity and `temperature` for temperature. 
    Currently, no method is provided to restore a flow variable by specifying the field name for the variable. 

    Note that all input parameters must be non-dimensionalized by the corresponding representative parameters. 
    PyPTS has no way to determine if the dimensions of these parameters are appropriate. 
    The user is solely responsible for the adequacy of the parameters specified. 

    """
    b_filename = filename.encode()
    if not read_ascii:
        raise NotImplementedError("read_binary is not implemented.")
    
    _pypts._pypts_flib.init_field_vtk(b_filename,
                                      ctypes.byref(ctypes.c_bool(read_ascii)))

    _pypts.logger.info("field data succesfully initialized from {0}".format(filename))

def init_field_afdet(filename:str):
    """
    Initializes the flow field. 

    The flow field is constructed from grid data stored in afdet solver format (.bin). 
    If the grid data also contains field data such as velocity and temperature fields, they are restored as well. 

    Parameters
    -----------
    filename (str) : filename in which the grid data is stored.

    Note
    ------
    Note that all input parameters must be non-dimensionalized by the corresponding representative parameters. 
    PyPTS has no way to determine if the dimensions of these parameters are appropriate. 
    The user is solely responsible for the adequacy of the parameters specified. 
    """
    b_filename = filename.encode()
    _pypts._pypts_flib.init_field_afdet(b_filename)

    _pypts.logger.info("field data succesfully initialized from {0}".format(filename))
