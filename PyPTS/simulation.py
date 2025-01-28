import numpy as np
import ctypes
from . import _pypts
import signal
import os

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



def initialize(nwrite:int, write_ascii:bool, traj_dir:str, 
               nback:int, backup_ascii:bool, backup_dir:str):
    """
    Initialize the simulation. 

    When initializing the simulation, the flow field, particle data, motion objects, and updaters 
    do not necessarily have to be initialized. 

    The particle data at each time are stored in two independent files: trajectory and backup. 

    The trajectory file (.vtk) is a Legacy VTK format file and is mainly used for visualization. 
    The position and state of a particle at each time are stored at each independent time step. 
    The user can control the time step interval at which the trajectory files are saved and the directory in which the files are saved. 
    The name of the saved file is `traj_path/trajectory__####.vtk`, 
    where `traj_path` is the path where the trajectory file is saved and #### is the time step where the data was recorded.

    The backup file (.pdata) contains all the data of the particles at each time step. 
    
    Parameters
    -----------
    nwrite (int) : Time step interval at which the particle trajectory visualization files (.vtk) will be output.
    write_ascii (bool) : Save trajectory files in ASCII.
    traj_dir (str) : Path of the directory where  trajectory files are stored. 
    nback (int)  : Time step interval at which the particle trajectory backup data (.pdata) will be output.
    backup_ascii (bool) : Save backupdata in ASCII.
    backup_dir (str) : Path of the directory where backup data ais stored. 
    """
    traj_path = os.path.join(traj_dir, "trajectory")
    backup_path = os.path.join(backup_dir, "particle")
    _pypts._pypts_flib.initialize_simulation(ctypes.byref(ctypes.c_int(nwrite)),
                                             ctypes.byref(ctypes.c_bool(write_ascii)),
                                             traj_path.encode(),
                                             ctypes.byref(ctypes.c_int(nback)),
                                             ctypes.byref(ctypes.c_bool(backup_ascii)),
                                             backup_path.encode()
                                             )
    
def run(nstart:int, nend:int):
    """
    Runs the simulation for the specified timesteps between `nstart` and `nend`. 

    At the time of execution, all the necessary data for the simulation (particle data, flow field data, updater, and motion) 
    must be prepared. 
    """
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    _pypts._pypts_flib.run_simulation(ctypes.byref(ctypes.c_int(nstart)),
                                      ctypes.byref(ctypes.c_int(nend)),
                                     )
    