import numpy as np
cimport numpy as np
cimport cython
from libcpp.map cimport map as cpp_map
# distutils: language = c++
np.import_array()


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef collapse_bins(np.ndarray bins, long binsize):
    cdef cpp_map[long, long] mapper
    cdef cpp_map[long, long].iterator it
    cdef long pos1, pos2
    cdef long pos
    cdef long size
    cdef long value
    cdef long[:,:] bin_view = bins.view(long)
    cdef long row_pos

    for row_pos in range(bin_view.size // 2):
        pos1, pos2 = bin_view[row_pos, :]
        for pos in range(pos1 + binsize, pos2, binsize):
            it = mapper.find(pos)
            if it == mapper.end():
                value = 1
            else:
                value = mapper[pos] + 1
            mapper[pos] = value

    return mapper
