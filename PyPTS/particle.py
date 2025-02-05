import numpy as np
import ctypes
from . import _pypts

class Particle(ctypes.Structure):
    """
    A struct representing a particle 
    
    It is compatible with `particle_t` in fortran `ParticleData.f90`

    Attributes:

    - pos      (3*double) : position
    - vel      (3*double) : velocity
    - f        (3*double) : force exerted on a particle
    - radius   (  double) : radius
    - ref_cell (     int) : cell id that a particle refers
    - state    (     int) : particle state (0: inactive, 1:active)
    
    """
    _fields_ = [
        ("pos", ctypes.c_double*3),
        ("vel", ctypes.c_double*3),
        ("f", ctypes.c_double*3),
        ("radius", ctypes.c_double),
        ("ref_cell", ctypes.c_int),
        ("state", ctypes.c_int),
        ("radius_min", ctypes.c_double),
        ("dead_time", ctypes.c_double)     
    ]


_pypts._pypts_flib.init_particle_data.argtypes = [ctypes.POINTER(ctypes.c_int),
                                                  np.ctypeslib.ndpointer(dtype=Particle)]
_pypts._pypts_flib.init_particle_data.restype = None


_pypts._pypts_flib.get_particle_data_size.argtypes = []
_pypts._pypts_flib.get_particle_data_size.restype = ctypes.c_int


_pypts._pypts_flib.get_particle_data.argtypes = [ctypes.POINTER(ctypes.c_int),
                                                  np.ctypeslib.ndpointer(dtype=Particle)]
_pypts._pypts_flib.get_particle_data.restype = None

_pypts._pypts_flib.restore_particle_data.argtypes = [ctypes.c_char_p]
_pypts._pypts_flib.restore_particle_data.restype  = None


class ParticleData:
    """
    particle data object.

    This is an interface of `mv_pdata` (`particle_data_t`). 

    Users can create data and pass it to fortran, or restore data from backup and copy it to Python.

    """
    def __init__(self, arg, copy=False) -> None:
        """
        Create a particle data. 

        
        A copy of particle data can be passed to Fortran. 

        Parameter
        ------------

        arg (int or str) : number of particles or filename of backup data.

        copy (bool)      : copy particle data from Fortran when load backup data.
        Note that  restored data is directly passed to Fortran and not passed to Python if `copy` is false.

        """
        if isinstance(arg, int):
            self.n_ = arg
            self.particles_ = np.empty(self.n_, dtype=Particle)
        elif isinstance(arg, str):
            _pypts._pypts_flib.restore_particle_data(arg.encode())
            _pypts.logger.info(f"Particle Data is restored from '{arg}'")
            if copy:
                self.get_particle_data()

        else:
            raise ValueError("Argument must be int or str.")

    @property
    def pos(self):
        """position of particles"""
        return self.particles_[:]["pos"]

    @pos.setter
    def pos(self, position:np.ndarray):
        self.particles_[:]["pos"] = position.reshape([self.n_,3]).astype(np.float64)


    @property
    def vel(self):
        """velocity of particles"""
        return self.particles_[:]["vel"]

    @vel.setter
    def vel(self, velocity:np.ndarray):
        self.particles_[:]["vel"] = velocity.reshape([self.n_,3]).astype(np.float64)

    @property
    def force(self):
        """forces exterted on particles"""
        return self.particles_[:]["f"]

    @force.setter
    def force(self, force:np.ndarray):
        self.particles_[:]["f"] = force.reshape([self.n_,3]).astype(np.float64)

    @property
    def radius(self):
        """radius of particles"""
        return self.particles_[:]["radius"]

    @radius.setter
    def radius(self, radius:np.ndarray):
        self.particles_[:]["radius"] = radius.reshape([self.n_]).astype(np.float64)

    @property
    def ref_cell(self):
        """cell ids that particles refer to"""
        return self.particles_[:]["ref_cell"]

    @ref_cell.setter
    def ref_cell(self, ref_cell:np.ndarray):
        self.particles_[:]["ref_cell"] = ref_cell.reshape([self.n_]).astype(np.int32)

    @property
    def state(self):
        """state of particles"""
        return self.particles_[:]["state"]

    @state.setter
    def state(self, state:np.ndarray):
        self.particles_[:]["state"] = state.reshape([self.n_]).astype(np.int32)

    @property
    def radius_min(self):
        """minimum radius of particles"""
        return self.particles_[:]["radius_min"]

    @radius_min.setter
    def radius_min(self, radius_min:np.ndarray):
        self.particles_[:]["radius_min"] = radius_min.reshape([self.n_]).astype(np.float64)

    @property
    def dead_time(self):
        """dead time of particles"""
        return self.particles_[:]["dead_time"]

    @dead_time.setter
    def dead_time(self, dead_time:np.ndarray):
        self.particles_[:]["dead_time"] = dead_time.reshape([self.n_]).astype(np.float64)


    def initialize_particle_data(self):
        """
        pass particle data to fortran. 

        """
        _pypts._pypts_flib.init_particle_data(ctypes.byref(ctypes.c_int(self.n_)),
                                              self.particles_)
        _pypts.logger.info("particle data passed to `init_particle_data`")

    def get_particle_data(self):
        """
        get current particle data from Fortran 
        """

        self.n_ = _pypts._pypts_flib.get_particle_data_size()
        self.particles_ = np.empty(self.n_, dtype=Particle)
        _pypts._pypts_flib.get_particle_data(ctypes.byref(ctypes.c_int(self.n_)), self.particles_)
        _pypts.logger.info("Get current particle data from Fortran")
        