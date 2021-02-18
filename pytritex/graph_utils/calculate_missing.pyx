import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
# distutils: language = c++

cdef extern from "math.h":
    float INFINITY


cdef np.ndarray[long, ndim=1, mode="c"] c_calculate_missing(
        np.ndarray[double, ndim=2, mode="c"] edges, np.ndarray[long, ndim=1, mode="c"] path):

    cdef long pos = 0
    cdef long opos = 0
    cdef double[:, :] edge_array = edges
    cdef long[:] path_array = path
    cdef vector[long] left_missing
    cdef vector[long] right_missing
    cdef vector[double] left_cost
    cdef vector[double] right_cost

    for pos in range(edges.shape[0]):
        if pos in path_array:
            continue
        for opos in range(1, path.shape[0] - 1):
            node = path[opos]
            if edges[pos, node] != np.nan:
                left_missing.pushback(pos)
                right_missing.pushback(opos)
                left_cost.pushback(edges[pos, node])
                right_cost.pushback(edges[node, path_array[opos + 1]])

