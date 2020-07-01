import numpy as np
cimport numpy as np
from libc.math cimport isnan
#distutils: language=c++
cdef extern from "math.h":
    float INFINITY
cimport cython
from cython.operator cimport dereference, postincrement
ctypedef pair[long, long] pll

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double _get_splicing_cost(double[:, :]& edges, cpp_map[long, pll]& pairs, long node) nogil:

    cdef:
        double splicing_cost = 0
        long previous
        long following
        double _cost

    previous = pairs[node].first
    following = pairs[node].second

    if previous > -1:
        # assert (edges[previous, node] > 0 and not np.isnan(edges[previous, node]))
        splicing_cost -= edges[previous, node]
    if following > -1:
        splicing_cost -= edges[node, following]

    if previous > -1 and following > -1:
        _cost = edges[previous, following]
        if not isnan(_cost) and _cost > 0:
            splicing_cost += _cost
        else:
            splicing_cost = INFINITY
    return splicing_cost


@cython.boundscheck(False)
@cython.wraparound(False)
cdef (double, long)  _get_insertion_cost(double[:, :]& edges,
                                         cpp_map[long, pll]& pairs,
                                         long header,
                                         long node) nogil:

    cdef:
        long opt_pos
        long other_node
        double cost, insertion_cost
        cpp_map[long, pll].iterator it = pairs.begin()
        double tmp_cost
        double first_cost, second_cost, third_cost
        bint first_valid, second_valid


    opt_pos = -1
    # Establish the first possible cost: inserting the node before every other node
    if node != header:
        cost = edges[header, node]
    else:
        cost = edges[node, pairs[node].second]

    if cost == 0 or isnan(cost):
        insertion_cost = INFINITY
    else:
        insertion_cost = cost
    # Now calculate the potential insertion cost for each possibility
    while it != pairs.end():
        other_node = dereference(it).first
        postincrement(it)
        if other_node == node or pairs[other_node].first == node or node == header:
            continue
        first_cost = edges[other_node, node]
        first_valid = first_cost != 0 and not isnan(first_cost)
        if pairs[other_node].second != -1: # If the node has a follower
            second_cost = edges[node, pairs[other_node].second]
            second_valid = not (second_cost == 0 or isnan(second_cost))
            third_cost = edges[other_node, pairs[other_node].second]
        else:
            second_cost = third_cost = 0
            second_valid = True

        if not (first_valid and second_valid):
            tmp_cost = INFINITY
        else:
            tmp_cost = first_cost + second_cost - third_cost

        if tmp_cost < insertion_cost:
            insertion_cost = tmp_cost
            opt_pos = other_node

    return insertion_cost, opt_pos


@cython.boundscheck(False)
@cython.wraparound(False)
def node_relocation(np.ndarray[DTYPE_f, ndim=2, mode="c"] edges,
                    np.ndarray[DTYPE_i, ndim=1, mode="c"] path,
                    double current_upper_bound,
                    long maxiter):

    # Edges is a SPARSE matrix where each non-found edge has a weight of 0
    # Path is a 1D numpy array with the positions.

    cdef:

        bint ever_changed
        bint changed
        long[:] path_view
        double[:, :] edge_view
        cpp_map[long, pll] pairs
        cpp_map[long, pll].iterator edge_it
        long position
        long header
        long pos
        long node
        double splicing_cost
        double insertion_cost
        long opt_pos
        long previous, following
        long new_following
        long iteration

    path_view = path[:]
    edge_view = edges[:]
    ever_changed = False
    changed = True
    pairs = cpp_map[long, pll]()
    edge_it = pairs.begin()
    # Create a map of positions.
    for position in range(path_view.shape[0]):
        node = path_view[position]
        pairs[node] = pll(-1, -1)
        if position > 0:
            pairs[node].first = path_view[position - 1]
        if position < path_view.shape[0] - 1:
            pairs[node].second = path[position + 1]

    header = path[0]
    iteration = 0
    while changed and iteration <= maxiter:
        changed = False
        iteration += 1
        for pos in range(path_view.shape[0]):
            # The *splicing* cost (given previous, node, next) is equal to:
            # weight[previous, next] - weight[node, next] - weight[previous, node]
            node = path_view[pos]
            splicing_cost = _get_splicing_cost(edges, pairs, node)
            insertion_cost, opt_pos = _get_insertion_cost(edges, pairs, header, node)
            if splicing_cost + insertion_cost < 0:
                # Moving the node costs less than keeping it in place. Exchange.
                ever_changed = True
                changed = True
                previous = pairs[node].first
                following = pairs[node].second
                if previous > -1:
                    pairs[previous].second = following
                if following > -1:
                    pairs[following].first = previous
                if opt_pos == -1:  # It's at the start
                    pairs[header].first = node
                    pairs[node] = pll(-1, header)
                    header = node
                else:
                    new_following = pairs[opt_pos].second
                    pairs[opt_pos].second = node
                    pairs[node] = pll(opt_pos, new_following)
                    pairs[new_following].first = node
                current_upper_bound += splicing_cost + insertion_cost  # This will result in a decrease

    if ever_changed:
        path[0] = header
        following = pairs[header].second
        for position in range(1, path_view.shape[0]):
            path[position] = following
            following = pairs[following].second

    return path, ever_changed, current_upper_bound
