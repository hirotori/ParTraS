import numpy as np
import ctypes
from . import _pypts

_pypts._pypts_flib.update_field_vtk.argtypes = [ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.c_char_p,
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_bool),
                                         ctypes.POINTER(ctypes.c_bool),
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_bool)]
_pypts._pypts_flib.update_field_vtk.restype = None

_pypts._pypts_flib.no_update_field.argtypes = None
_pypts._pypts_flib.no_update_field.restype = None

# TODO: updaterのベースクラスを作って, それを継承する形にする. 
#       今後ほかのファイルに対するupaderが必要になるので. 

def no_update():
    """
    no update event occured.
    """
    _pypts._pypts_flib.no_update_field()
    _pypts.logger.info("updater disabled")


class vtkUpdater:
    def __init__(self, dt_f:float, dt_p:float, basename:str, ndigit:int, field_only:bool, verts_only:bool, 
                 interval:int, ascii:bool) -> None:
        """
        create interface object corresponding with the field updater instance in Fortran.

        Assume that any flow field file at time step `ncyc_flow` is named as `basename+"_"+str(ncyc_flow)+".vtk"`. 
        if ndigit > 0, find filename `basename+"_"+str(ncyc_flow).zfill(ndigit)+".vtk"`

        Parameters
        -----------
        dt_f (float) : time interval of flow field
        dt_p (float) : time interval of particle simulation
        basename (str) : filename of flow field without extension
        ndigit (int) : number of digit in filename. if digit = 0, no zero-padding exists.
        field_only (bool) : force to update only field variable (velocity). If False, update all members of flow field.
        verts_only (bool) : force to update only vertices of cells. This flag is valid only if field_only is `False`.
        interval (int) : output interval of flow field
        ascii (bool) : read file in ASCII format. if False, read in binary
        """
        self.dt_f = dt_f
        self.dt_p = dt_p
        self.basename = basename
        self.ndigit = ndigit
        self.field_only = field_only
        self.verts_only = verts_only
        self.interval = interval
        self.ascii = ascii

        self.enable()

    def disable(self):
        """
        disable updater. 
        importers associated with updater are deallocated. 
        """
        no_update()

    def enable(self):
        """
        enable updater

        Example: 
        
        ```Python

        # automatically enabled
        updater = vtkUpdater(dt_f, dt_p, basename, pad, field_only, interval, ascii)
        # ...
        updater.disable()
        # ...
        # reload
        updater.enable()
        
        ```
        """
        _pypts._pypts_flib.update_field_vtk(ctypes.byref(ctypes.c_double(self.dt_f)),
                                            ctypes.byref(ctypes.c_double(self.dt_p)),
                                            self.basename.encode(),
                                            ctypes.byref(ctypes.c_int(self.ndigit)),
                                            ctypes.byref(ctypes.c_bool(self.field_only)),
                                            ctypes.byref(ctypes.c_bool(self.verts_only)),
                                            ctypes.byref(ctypes.c_int(self.interval)),
                                            ctypes.byref(ctypes.c_bool(self.ascii)))