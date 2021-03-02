cimport numpy as np
cimport cython
from libcpp.map cimport map as cpp_map

cpdef cpp_map[long,long]  collapse_bins(np.ndarray[dtype=long, ndim=2] bins, long binsize)  except +