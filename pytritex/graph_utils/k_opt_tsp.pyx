import numpy as np
cimport numpy as np

ctypedef np.int_t DTYPE_int
cdef extern from "math.h":
    float INFINITY


cdef double c_route_cost(np.ndarray[double, ndim=2, mode="c"] graph, np.ndarray[long, ndim=1, mode="c"] path):
    cdef double cost = 0
    cdef long index = 0
    cdef double[:, :] graph_array = graph
    cdef long[:] path_array = path
    cdef long shape = path_array.shape[0]
    cdef double temp_cost
    cost = graph_array[path_array[shape - 1], path_array[0]]

    for index in range(shape - 1):
        temp_cost = graph_array[path_array[index]][path_array[index + 1]]
        if temp_cost == 0:
            return INFINITY
        else:
            cost += graph_array[path_array[index]][path_array[index + 1]]

    return cost


cpdef float route_cost(np.ndarray graph, np.ndarray path):

    cdef np.ndarray[double, ndim=2, mode="c"] graph_array = graph.astype(float)
    cdef np.ndarray[long, ndim=1, mode="c"] path_array = path
    return c_route_cost(graph_array, path_array)


cdef np.ndarray[long, ndim=1, mode="c"] _swap_2opt(np.ndarray[long, ndim=1, mode="c"] route, long i, long k):
    """ Swapping the route."""

    cdef long[:] route_view = route
    cdef np.ndarray[DTYPE_int, ndim=1, mode="c"] new_route = np.empty(route_view.shape[0], dtype=np.int)
    new_route = np.concatenate([route[:i], route[i:k + 1][::-1], route[k + 1:]])
    return new_route


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
    cdef np.ndarray[double, ndim=2, mode="c"] graph_array = graph.astype(float)
    cdef np.ndarray[long, ndim=1, mode="c"] route_array = route
    cdef np.ndarray[long, ndim=1, mode="c"] best_found_route
    cdef np.ndarray[long, ndim=1, mode="c"] new_route
    cpdef long[:] best_found_route_view

    best_found_route = route_array[:]
    best_found_route_view = best_found_route
    max_index = best_found_route_view.shape[0] - 1
    best_cost = c_route_cost(graph_array, best_found_route)
    while improved == 1:
        improved = 0
        for index in range(1, max_index):
            for kindex in range(index + 1, max_index):
                new_route = _swap_2opt(best_found_route, index, kindex)
                cost = c_route_cost(graph_array, new_route)
                if cost < best_cost:
                    best_cost = cost
                    best_found_route = new_route
                    improved = 1
            if improved:
                break
    return best_found_route
