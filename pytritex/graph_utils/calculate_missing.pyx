import numpy as np
cimport numpy as np
cimport cython
from libcpp.vector cimport vector
# distutils: language=c++

cdef extern from "math.h":
    float INFINITY
    bint isnan(double x)


@cython.boundscheck(False)
cdef bint in_array(long val, vector[long] arr):
    cdef int pos
    cdef long size = arr.size()
    for pos in range(size):
        if arr[pos] == val:
            return True
    return False


@cython.boundscheck(False)
cdef np.ndarray[long, ndim=1, mode="c"] c_calculate_missing(
        np.ndarray[double, ndim=2, mode="c"] edges, np.ndarray[long, ndim=1, mode="c"] path):

    cdef long pos = 0
    cdef long opos = 0
    cdef double[:, :] edge_array = edges
    cdef vector[long] path_array = path
    cdef vector[long] left_missing
    cdef vector[long] right_missing
    cdef vector[double] left_cost
    cdef vector[double] right_cost
    cdef long node
    cdef double cost

    for pos in range(edges.shape[0]):
        if in_array(pos, path_array):
            continue
        for opos in range(1, path.shape[0] - 1):
            node = path_array[opos]
            if not isnan(edges[pos, node]):
                left_missing.push_back(pos)
                right_missing.push_back(opos)
                cost = edge_array[pos, node]
                left_cost.push_back(cost)
                cost = edge_array[node, path_array[opos + 1]]
                right_cost.push_back(cost)

