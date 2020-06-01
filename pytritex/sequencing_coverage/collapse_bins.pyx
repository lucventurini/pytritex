import numpy as np
cimport numpy as np
cimport cython
from libcpp.map cimport map as cpp_map
# distutils: language=c++
np.import_array()


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef collapse_bins(np.ndarray bins, long binsize):
    cdef cpp_map[long, long] mapper
    cdef cpp_map[long, long].iterator it
    cdef long pos1, pos2
    cdef long pos, second_pos
    cdef int size
    cdef int value
    cdef np.ndarray[np.int_t, ndim=2, mode="c"] bin_view = bins
    cdef int row_pos
    cdef long binsize_c = binsize

    for row_pos in range(0, bin_view.shape[0], 1):
        pos1 = bin_view[row_pos, 0] + binsize_c
        pos2 = bin_view[row_pos, 1]
        while pos < pos2:
            it = mapper.find(pos)
            if it == mapper.end():
                value = 1
            else:
                value = mapper[pos] + 1
            mapper[pos] = value
            pos += binsize_c

    return mapper
