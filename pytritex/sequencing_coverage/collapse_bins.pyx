import numpy as np
cimport numpy as np
cimport cython
from libcpp.map cimport map as cpp_map
# distutils: language=c++
np.import_array()


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef collapse_bins(np.ndarray bins, long binsize, long num_threads):
    cdef cpp_map[long, long] mapper
    cdef cpp_map[long, long].iterator it
    cdef long pos1, pos2
    cdef int size
    cdef int value
    cdef long steps, step
    cdef long[:, :] bin_view = bins
    cdef long row_pos
    cdef long binsize_c = binsize
    cdef long nrows = bin_view.shape[0]

    for row_pos in range(0, nrows, 1):
        pos1 = bin_view[row_pos, 0]
        pos2 = bin_view[row_pos, 1]
        for step in range(1, (pos2 - pos1) // binsize):
            if pos1 + step * binsize_c >= pos2:
                break
            it = mapper.find(pos1 + step * binsize_c)
            if it == mapper.end():
                value = 1
            else:
                value = mapper[pos1 + step * binsize_c] + 1
            mapper[pos1 + step * binsize_c] = value

    return mapper
