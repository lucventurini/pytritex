import numpy as np
np.seterr(all="raise")

# TRAVELING SALESMAN PROBLEM
# Insertion Algorithms
# (Rosenkrantz, Stearns, Lewis, 1974)
#
# An insertion procedure takes a sub-tour on k nodes at iteration k and determines which of the remaining n-k nodes shall be inserted to the sub-tour next (the selection step) and where (between which two nodes) it should be inserted (the insertion step).
#
# Nearest Insertion
#
#     Step 1. Start with a sub-graph consisting of node i only.
#     Step 2. Find node r such that cir is minimal and form sub-tour i-r-i.
#     Step 3. (Selection step) Given a sub-tour, find node r not in the sub-tour closest to any node j in the sub-tour; i.e. with minimal crj
#     Step 4. (Insertion step) Find the arc (i, j) in the sub-tour which minimizes cir + crj - cij . Insert r between i and j.
#     Step 5. If all the nodes are added to the tour, stop. Else go to step 3


def insert_nodes(edges: np.ndarray, path: np.ndarray):

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

    edges = edges[:]
    visited = dict((_, _ in path) for _ in range(edges.shape[0]))
    header = path[0]
    dpath = dict()
    for pos, el in enumerate(path):
        link = [-1, -1]
        if pos > 0:
            link[0] = path[pos - 1]
        if pos < path.shape[0] - 1:
            link[1] = path[pos + 1]
        dpath[el] = link
        visited[el] = True

    for vertex, value in visited.items():
        if value is True:
            continue
        visited[vertex] = True
        current_position = header
        opt_position = -1
        opt_cost = edges[header, vertex]
        if np.isnan(opt_cost):
            opt_cost = float("inf")
        while current_position != -1:
            cost = edges[current_position, vertex]
            if not np.isnan(cost):
                following = dpath[current_position][1]
                if following != -1:
                    add_cost = edges[vertex, following]
                    if np.isnan(add_cost):
                        cost = float("inf")
                    else:
                        remove_cost = edges[current_position, following]
                        cost = cost + add_cost - remove_cost
            else:
                cost = float("inf")
            if cost < opt_cost:
                opt_cost = cost
                opt_position = current_position
            current_position = dpath[current_position][1]
            continue

        if opt_cost < float("inf"):  # Only insert nodes for which there is a link.
            if opt_position == -1:
                dpath[vertex] = [-1, header]
                dpath[header][0] = vertex
                header = vertex
            else:
                new_following = dpath[opt_position][1]
                dpath[vertex] = [opt_position, new_following]
                dpath[opt_position][1] = vertex
                if new_following != -1:
                    dpath[new_following][0] = vertex

    path = [header]
    current_iter = 0
    current_position = header
    while current_iter < len(dpath):
        following = dpath[current_position][1]
        if following == -1:
            break
        path.append(following)
        current_position = following
    assert len(path) == len(dpath)
    path = np.array(path)
    assert np.unique(path).shape[0] == path.shape[0]
    return path
