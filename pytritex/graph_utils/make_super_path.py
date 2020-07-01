import pandas as pd
import numpy as np
import networkit as nk
from sys import float_info
import scipy.sparse.csgraph
import scipy.sparse
from .k_opt_tsp import tsp_2_opt
from .insert_nodes import insert_nodes
from .node_relocation import node_relocation
import scipy.stats
import scipy.special
scipy.special.seterr(all="raise")


# This is necessary for networkit, as currently the distances can be "infinite"
# See https://github.com/networkit/networkit/issues/541
max_double = float_info.max


def convert_array(from_keys, to_keys, original):
    sort_idx = np.argsort(from_keys)
    idx = np.searchsorted(from_keys, original, sorter=sort_idx)
    out = to_keys[sort_idx][idx]
    return out


def _assign_ranks(edges, path, cidx, backbone):
    """This method will assign to each node its *rank*, ie *its distance from the backbone*
    as measured by the number of edges."""

    # Create a vertical table. Each row is a node: node, bin, rank, backbone
    ranks = np.repeat(np.arange(cidx.shape[1], dtype=np.float), 4).reshape(cidx.shape[1], 4)
    ranks[:, [1, 2]] = np.nan
    # Assign the rank (0) and the bin
    ranks[path, 1] = np.arange(path.shape[0])
    ranks[path, 2] = 0
    rank = 0
    ranks[:, 3] = False
    ranks[backbone, 3] = True

    foundany = True
    while np.isnan(ranks).any():
        rank += 1
        foundany = False
        curr_nan = np.where(np.isnan(ranks[:, 1]))[0]
        curr_not_na = np.where(~np.isnan(ranks[:, 1]))[0]
        for node in curr_nan:
            edge = edges[node]
            cols = np.intersect1d(np.where(~np.isnan(edge)), curr_not_na, assume_unique=True)
            if cols.shape[0] > 0:
                best_score = edge[cols].max()
                best_node = cols[np.where(edge[cols] == best_score)][0]
                ranks[node, [1, 2]] = ranks[best_node, 1], rank
                foundany = True
        if not foundany:
            # Infinite loop otherwise! Some nodes might be disconnected
            break

    nulls = np.isnan(ranks[:, 2])
    ranks[nulls, 1] = -1
    ranks[nulls, 2] = max(2, rank + 1)
    ranks[nulls, 3] = False
    return ranks


def _get_path(matrix, cidx, start=None, end=None, ncores=1):
    """Initial path. This method will:
    - calculate the minimum spanning tree from the edges
    - calculate the longest simple path (or the shortest path between start and end if given)
    - return: the path, the minimum spanning tree, and the total edge cost of the MST.
    """
    mst = scipy.sparse.csgraph.minimum_spanning_tree(matrix)

    # This is our minimum cost after calculating the MST
    # the sparse matrix from nk will already contain the weights
    # As the matrix is symmetrical, it's the sum of the weights divided by 2.
    mst_lower_bound = mst.toarray().sum() / 2
    # Now back to NK
    nk.setNumberOfThreads(ncores)
    mgraph = nk.Graph(mst.shape[0], weighted=True, directed=False)
    marr = matrix.toarray()
    for x, y in zip(mst.nonzero()[0], mst.nonzero()[1]):
        if marr[x, y] == 0:
            continue
        mgraph.addEdge(x, y, 1)
    # Now that we have the new minimum-spanning-tree graph we can get the diameter or the distance
    nodes = np.fromiter(mgraph.iterNodes(), np.int)
    if start is None or end is None:
        # This runs the Dijkstra's algorithm of BFS iteratively on each node.
        apsp = nk.distance.APSP(mgraph).run()
        distances = np.array(apsp.getDistances())
        # Set unavailable distances to -1
        distances[(distances == max_double) | np.isnan(distances)] = -1
        best = np.where(distances == distances.max())
        start, end = nodes[best[0][0]], nodes[best[1][0]]
    else:
        # Get the local coordinates
        start = cidx[:, np.where(cidx[0, :] == start)[0]][1]
        end = cidx[:, np.where(cidx[0, :] == end)[0]][1]

    # Reconstruct the backbone path.
    path = np.array(nk.distance.BFS(mgraph, start, target=end).run().getPath(end), dtype=np.int)

    return path, mst, mst_lower_bound


def _local_improvement(edges, path, mst, current_upper_bound, mst_lower_bound, ncores, maxiter):
    if current_upper_bound > mst_lower_bound:
        changed = True
        rounds = 0
        while changed:
            internal_rounds = 0
            while changed:
                prev = path[:]
                path = tsp_2_opt(mst.toarray(), path, ncores)
                if all(path == prev):
                    changed = False
                else:
                    changed = True
                # The current upper bound is now the new weight cost
                current_upper_bound = sum(
                    edges[path[index], path[index + 1]] for index in np.arange(path.shape[0] - 1,
                                                                                    dtype=np.int))
                internal_rounds += 1
                if internal_rounds >= maxiter:
                    break
            path, changed, current_upper_bound = node_relocation(edges, path, current_upper_bound)
            rounds += 1
            if rounds >= maxiter:
                break

    return path, current_upper_bound


# TODO: remove all unnecessary pandas code from here. It only slows down the function.
def make_super_path(origin_edge_list, cms, start=None, end=None, maxiter=100, verbose=True, ncores=1):
    """Order scaffolds using Hi-C links for one chromosome.
    This procedure is based on Y. Wu et al. 2008,  doi:10.1371/journal.pgen.1000212

    """

    # Get backbone from minimum spanning tree (MST)
    origin_edge_list = origin_edge_list[["cluster1", "cluster2", "weight"]].values
    origin_edge_list = origin_edge_list[origin_edge_list[:, 0] != origin_edge_list[:, 1]][:]
    # Now make it symmetrical
    origin_edge_list = np.vstack(
        [origin_edge_list,
         np.vstack([origin_edge_list[:, 1], origin_edge_list[:, 0], origin_edge_list[:, 2]]).T])

    orig_coords, weights = origin_edge_list[:, [0, 1]], origin_edge_list[:, 2]
    # Transpose the ids to a new system of coordinates.
    cidx = np.unique(orig_coords.flatten())
    # Column *0*: original scaffold indices
    # Column *1*: new partial scaffold indices (from 0 to shape[0] - 1, inclusive)
    cidx = np.vstack([cidx, np.arange(cidx.shape[0])])
    coords = convert_array(cidx[0, :], cidx[1,:], orig_coords)
    # Assign the new index
    try:
        shape = [int(coords.max() + 1)] * 2
    except Exception:
        np.save("/tmp/cidx.arr", cidx)
        np.save("/tmp/orig_coords.arr", orig_coords)
        np.save("/tmp/coords.arr", coords)
        pd.to_pickle(cms, "/tmp/cms.pkl")
        pd.to_pickle(origin_edge_list, "/tmp/oel.pkl")
        raise

    assert weights.min() > 0
    matrix = scipy.sparse.coo_matrix((weights, (coords[:, 0], coords[:, 1])),
                                     dtype=weights.dtype, shape=shape)
    edges = matrix.toarray()
    # Is the matrix symmetrical? If it is not, we have a problem.
    assert (edges - edges.T == 0).all(), (edges.tolist())
    edges[edges == 0] = np.nan
    try:
        path, mst, mst_lower_bound = _get_path(matrix, cidx, start=start, end=end, ncores=ncores)
    except AssertionError:
        raise AssertionError(orig_coords)
    for pos in range(len(path) - 1):
        assert not np.isnan(edges[path[pos], path[pos + 1]])
        assert not np.isnan(edges[path[pos + 1], path[pos]])

    # These are the nodes that constitute the backbone. Keep them aside for now.
    backbone = path[:]
    assert np.in1d(path, cidx).all()
    _prev = path[:]
    path = insert_nodes(edges, path)
    # Now get the current upper bound. This is the sum of the current path.
    current_upper_bound = sum(edges[path[index], path[index + 1]]
                              for index in np.arange(path.shape[0] - 1, dtype=np.int))
    cost_after_initialization = current_upper_bound
    # Run the improvement procedure ONLY IF the CUB is greater than the MST cost.
    path, current_upper_bound = _local_improvement(edges, path, mst,
                                                   current_upper_bound, mst_lower_bound, ncores,
                                                   maxiter)
    # We will presume that we have to SKIP the block_optimise procedure
    # Assign to each node its rank, ie the distance from the backbone.
    ranks = _assign_ranks(edges, path, cidx, backbone)
    ranks = np.hstack(
        [convert_array(cidx[1, :], cidx[0, :], ranks[:, 0]).reshape(ranks.shape[0], 1),
         ranks[:, [1, 2, 3]]])

    try:
        cms = cms.reindex(ranks[:, 0])["cM"].values
        index = np.where(~np.isnan(cms))[0]
        if index.shape[0] > 1:
            _cms = cms[index]
            _bins = ranks[index, 1].reshape(index.shape[0])
            assert _cms.shape[0] == _bins.shape[0], (_cms, _bins)
            if np.unique(_cms).shape[0] < 2:
                flip = np.nan
            else:
                flip = scipy.stats.pearsonr(_bins, _cms)[0]
        else:
            flip = np.nan
    except KeyError as exc:
        print(cms.index)
        print(ranks[:, 0])
        print(cidx)
        raise KeyError(exc)

    if not np.isnan(flip) and flip < 0:
        # Flip the order
        ranks[:, 1] = ranks[:, 1].max() - ranks[:, 1] + 1
        # df.eval("bin = bin.max() - bin", inplace=True)

    return ranks
