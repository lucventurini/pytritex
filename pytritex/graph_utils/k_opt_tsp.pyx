import numpy as np
cimport numpy as np
cimport cython
np.import_array()

ctypedef np.int_t DTYPE_int
cdef extern from "math.h":
    float INFINITY
from cpython cimport array
import array


@cython.wraparound(False)
@cython.boundscheck(False)
cdef double c_route_cost(double[:, :] graph, long[:] path) nogil:
    cdef:
        double cost
        long shape = path.shape[0]
        double temp_cost
        long index = shape - 1
        long second = 0

    cost = graph[path[index], path[0]]
    for index in range(shape - 1):
        second = index + 1
        temp_cost = graph[path[index]][path[second]]
        if temp_cost == 0:
            return INFINITY
        else:
            cost += temp_cost

    return cost


cpdef float route_cost(np.ndarray graph, np.ndarray path):

    cdef np.ndarray[double, ndim=2, mode="c"] graph_array = graph.astype(float)
    cdef np.ndarray[long, ndim=1, mode="c"] path_array = path
    return c_route_cost(graph_array, path_array)

@cython.wraparound(False)
@cython.boundscheck(False)
cdef long[:] _swap(long[:] route_array, long index, long kindex):
    cdef:
        long[:] new_route
        long internal_index
        long mirror
        long[:] to_swap

    new_route = route_array.copy()
    to_swap = route_array[index:kindex + 1]
    cdef long swap_length = kindex - index + 1
    for internal_index in range(0, swap_length, 1):
        mirror = swap_length - internal_index - 1
        new_route[index + internal_index] = to_swap[mirror]
    return new_route


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray tsp_2_opt(np.ndarray graph, np.ndarray route):
    """
    Approximate the optimal path of travelling salesman according to 2-opt algorithm
    Args:
        graph: 2d numpy array as graph
        route: list of nodes

    Returns:
        optimal path according to 2-opt algorithm

    Examples:
        >>> import numpy as np
        >>> graph = np.array([[  0, 300, 250, 190, 230],
        >>>                   [300,   0, 230, 330, 150],
        >>>                   [250, 230,   0, 240, 120],
        >>>                   [190, 330, 240,   0, 220],
        >>>                   [230, 150, 120, 220,   0]])
        >>> tsp_2_opt(graph)
    """
    cdef bint improved = 1
    cdef double best_cost
    cdef double cost
    cdef long index
    cdef long kindex
    cdef long max_index
    cdef long route_shape
    cdef double[:, :] graph_array = graph.astype(float)
    cdef long[:] route_array = route
    cdef long[:] best_found_route = route_array[:]
    cdef long[:] new_route
    cdef long swap_length
    cdef np.ndarray[DTYPE_int, ndim=1] final_route

    max_index = best_found_route.shape[0] - 1
    best_cost = c_route_cost(graph_array, best_found_route)
    while improved == 1:
        improved = 0
        for index in range(1, max_index):
            for kindex in range(index + 1, max_index):
                # Swap internally between index and kindex
                new_route = _swap(route_array, index, kindex)
                cost = c_route_cost(graph_array, new_route)
                if cost < best_cost:
                    best_cost = cost
                    best_found_route = new_route
                    improved = 1
            if improved:
                break

    final_route = np.array(best_found_route[:])
    return final_route
