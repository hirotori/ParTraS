import os
import ctypes
from numpy import ctypeslib
import numpy as np

from logging import getLogger, StreamHandler, INFO, Formatter
logger = getLogger("PyPTS")
logger.setLevel(INFO)
hdlr = StreamHandler()
format = Formatter('%(asctime)s:(%(name)s)::%(levelname)s:: %(message)s')
hdlr.setFormatter(format)
hdlr.setLevel(INFO)
logger.addHandler(hdlr)

_libpath = os.path.dirname(os.path.abspath(__file__))
_pypts_flib = ctypeslib.load_library("libpypts.dylib", _libpath)

