import numpy as np
import ctypes
from . import _pypts
import signal

# init
_pypts._pypts_flib.initialize_simulation.argtypes = [ctypes.POINTER(ctypes.c_int),
                                                     ctypes.POINTER(ctypes.c_bool),
                                                     ctypes.c_char_p,
                                                     ctypes.POINTER(ctypes.c_int),
                                                     ctypes.POINTER(ctypes.c_bool),
                                                     ctypes.c_char_p
                                                     ]
_pypts._pypts_flib.initialize_simulation.restype = None


# run
_pypts._pypts_flib.run_simulation.argtypes = [ctypes.POINTER(ctypes.c_int),
                                              ctypes.POINTER(ctypes.c_int),
                                             ]
_pypts._pypts_flib.run_simulation.restype = None


# settings
_pypts._pypts_flib.set_writeout_settings.argtypes = [ctypes.POINTER(ctypes.c_int),
                                                     ctypes.POINTER(ctypes.c_bool),
                                                     ctypes.c_char_p
                                                     ]
_pypts._pypts_flib.set_writeout_settings.restype = None


# settings
_pypts._pypts_flib.set_dump_settings.argtypes = [ctypes.POINTER(ctypes.c_int),
                                                 ctypes.POINTER(ctypes.c_bool),
                                                 ctypes.c_char_p
                                                ]
_pypts._pypts_flib.set_dump_settings.restype = None



def initialize(nwrite:int, write_ascii:bool, traj_path:str, 
               nback:int, backup_ascii:bool, backup_path:str):
    _pypts._pypts_flib.initialize_simulation(ctypes.byref(ctypes.c_int(nwrite)),
                                             ctypes.byref(ctypes.c_bool(write_ascii)),
                                             traj_path.encode(),
                                             ctypes.byref(ctypes.c_int(nback)),
                                             ctypes.byref(ctypes.c_bool(backup_ascii)),
                                             backup_path.encode()
                                             )
    
def run(nstart:int, nend:int):
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    _pypts._pypts_flib.run_simulation(ctypes.byref(ctypes.c_int(nstart)),
                                      ctypes.byref(ctypes.c_int(nend)),
                                     )
    