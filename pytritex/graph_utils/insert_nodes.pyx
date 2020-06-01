import numpy as np
cimport numpy as np
np.seterr(all="raise")
np.import_array()
from libcpp.map cimport map as cpp_map
from libcpp.pair cimport pair
from libc.math cimport isnan
# distutils: language=c++
cdef extern from "math.h":
    float INFINITY
ctypedef pair[long, long] pll
cimport cython
from cython.operator cimport dereference, postincrement


@cython.boundscheck(False)
cdef (cpp_map[long, pll], long) c_insert_nodes(cpp_map[long, bint] visited, cpp_map[long, pll] dpath,
                                     double[:, :] edges, long header):

    cdef long current_position, opt_pos, new_following
    cdef long vertex
    cdef bint visit_check
    cdef double opt_cost, cost, remove_cost, add_cost
    cdef cpp_map[long, bint].iterator it = visited.begin()

    while (it != visited.end()):
        vertex = dereference(it).first
        visit_check = dereference(it).second
        postincrement(it)
        if visit_check == 1:
            continue
        visited[vertex] = 1
        current_position = header
        opt_position = -1
        opt_cost = edges[header, vertex]
        if isnan(opt_cost):
            opt_cost = INFINITY

        while current_position != -1:
            cost = edges[current_position, vertex]
            if isnan(cost):
                cost = INFINITY
            else:
                following = dpath[current_position].second
                if following != -1:
                    add_cost = edges[vertex, following]
                    if isnan(add_cost):
                        cost = INFINITY
                    else:
                        remove_cost = edges[current_position, following]
                        cost = cost + add_cost - remove_cost

            if cost < opt_cost:
                opt_cost = cost
                opt_position = current_position
            current_position = dpath[current_position].second
            continue

        if opt_cost < INFINITY:
            if opt_position == -1:
                dpath[vertex].first = -1
                dpath[vertex].second = header
                dpath[header].first = vertex
                header = vertex
            else:
                new_following = dpath[opt_position].second
                dpath[vertex].first = opt_position
                dpath[vertex].second = new_following
                dpath[opt_position].second = vertex
                if new_following != -1:
                    dpath[new_following].first = vertex

    return dpath, header


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[DTYPE_i, ndim=1] insert_nodes(np.ndarray[DTYPE_f, ndim=2, mode="c"] edges,
                                               np.ndarray[DTYPE_i, ndim=1, mode="c"] path):

    """Implementation of the TSP insertion algorithm (Rosenkrantz, Stearns, Lewis, 1974)
    to find the location of missing nodes from the longest path (ie the backbone).

    Given a path P, which is a *sub-tour* of the graph:
    1- Given a sub-tour, find node r *not* in the sub-tour closest to any node j in the sub-tour; i.e. with minimal crj
    2- Find the arc (i, j) in the sub-tour which minimizes cir + crj - cij . Insert r between i and j.
    3- If all the nodes are added to the tour, stop. Else go to step 3

    :param edges: a numpy array (table) of shape (N, N), with N being the number of vertices.
                  Each filled position in the table corresponds to an edge.
    :type edges: np.ndarray
    :param path: 1-dim numpy array with the correct order of edges.
    :type path: np.ndarray
    :returns: new edges and path.
    """

    # Edge list is a np.ndarray of shape [x, 3] with columns corresponding to
    # node1, node2, weight
    # path is an ordered list of nodes,

    cdef double[:, :] edges_view = edges
    cdef long[:] path_view = path.copy()
    cdef cpp_map[long, bint] visited
    cdef cpp_map[long, pll] dpath
    cdef pll link
    cdef long visitor
    cdef long header, el
    cdef long current_iter, following, current_position

    for visitor in range(edges_view.shape[0]):
        visited[visitor] = 0

    header = path_view[0]
    for visitor in range(path_view.shape[0]):
        el = path_view[visitor]
        link.first, link.second = -1, -1
        if visitor > 0:
            link.first = path_view[visitor - 1]
        if visitor < path_view.shape[0] - 1:
            link.second = path_view[visitor + 1]
        dpath[el] = link
        visited[el] = 1
    dpath, header = c_insert_nodes(visited, dpath, edges, header)

    path = np.empty(dpath.size(), dtype=np.int)
    current_iter = 0
    path[current_iter] = header
    current_position = header
    while current_iter < dpath.size():
        current_iter += 1
        following = dpath[current_position].second
        if following == -1:
            break
        path[current_iter] = following
        current_position = following

    return path
