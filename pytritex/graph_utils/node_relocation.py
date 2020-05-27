import numpy as np
np.seterr(all="raise")


def _get_splicing_cost(edges, pairs, node):
    splicing_cost = 0
    previous, following = pairs[node]
    if previous > -1:
        assert (edges[previous, node] > 0 and not np.isnan(edges[previous, node]))
        splicing_cost -= edges[previous, node]
    if following > -1:
        if not (edges[following, node] > 0 and not np.isnan(edges[following, node])):
            # Something has gone wrong here. The path should always have defined edges!
            print(pairs)
            print(node, following)
            print(edges)
            raise AssertionError

        splicing_cost -= edges[node, following]

    if previous > -1 and following > -1:
        _cost = edges[previous, following]
        if not np.isnan(_cost) and _cost > 0:
            splicing_cost += _cost
        else:
            splicing_cost = float("inf")
    return splicing_cost


def _get_insertion_cost(edges, pairs, header, node, position):
    opt_pos = -1
    # Establish the first possible cost: inserting the node before every other node
    if node != header:
        cost = edges[header, node]
    else:
        cost = edges[node, pairs[node][1]]
    if cost == 0 or np.isnan(cost):
        insertion_cost = float("inf")
    else:
        insertion_cost = cost
    # Now calculate the potential insertion cost for each possibility
    for other_node in pairs:
        if node == other_node or node == pairs[node][0] or node == header:
            # Skip self calculation or calculation with the previous index
            continue
        first_cost = edges[other_node, node]
        first_valid = not (first_cost == 0 or np.isnan(first_cost))
        if pairs[other_node][1] != -1:  # If the node has a follower
            second_cost = edges[node, pairs[other_node][1]]
            second_valid = not (second_cost == 0 or np.isnan(second_cost))
            third_cost = edges[other_node, pairs[other_node][1]]
            assert not (third_cost == 0 or np.isnan(third_cost))
        else:
            second_cost = third_cost = 0
            second_valid = True
        if not (first_valid and second_valid):
            tmp_cost = float("inf")
        else:
            tmp_cost = first_cost + second_cost - third_cost

        if tmp_cost < insertion_cost:
            insertion_cost = tmp_cost
            opt_pos = position

    return insertion_cost, opt_pos


def node_relocation(edges, path: np.array, current_upper_bound):

    # Edges is a SPARSE matrix where each non-found edge has a weight of 0
    # Path is a 1D numpy array with the positions.

    ever_changed = False
    changed = True
    pairs = dict()
    # Create a map of positions.
    for position in np.arange(path.shape[0], dtype=np.int):
        node = path[position]
        pairs[node] = [-1, -1]
        if position > 0:
            pairs[node][0] = path[position - 1]
        if position < path.shape[0] - 1:
            pairs[node][1] = path[position + 1]

    header = path[0]
    while changed:
        changed = False
        for pos in np.arange(path.shape[0]):
            # The *splicing* cost (given previous, node, next) is equal to:
            # weight[previous, next] - weight[node, next] - weight[previous, node]
            node = path[pos]
            assert node in pairs, (node, pairs)
            splicing_cost = _get_splicing_cost(edges, pairs, node)
            insertion_cost, opt_pos = _get_insertion_cost(edges, pairs, header, node, pos)
            if splicing_cost + insertion_cost < 0:
                # Moving the node costs less than keeping it in place. Exchange.
                ever_changed = changed = True
                previous, following = pairs[node]
                if previous > -1:
                    pairs[previous][1] = following
                if following > -1:
                    pairs[following][0] = previous
                if opt_pos == -1:  # It's at the start
                    pairs[header][0] = node
                    pairs[node] = [-1, header]
                    header = node
                else:
                    new_following = pairs[opt_pos][1]
                    pairs[opt_pos][1] = node
                    pairs[node] = [opt_pos, new_following]
                    pairs[new_following][0] = node
                current_upper_bound += splicing_cost + insertion_cost  # This will result in a decrease

    if ever_changed:
        path[0] = header
        following = pairs[header][1]
        for position in np.arange(1, path.shape[0]):
            path[position] = following
            following = pairs[following][1]

    return path, ever_changed, current_upper_bound
