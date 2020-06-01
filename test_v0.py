import joblib
import pandas as pd
import itertools
from pytritex.scaffold_10x.__init__ import scaffold_10x
from pytritex.utils import n50
import numexpr as ne
import itertools
import multiprocessing as mp
import _pickle
from collections import namedtuple
import dask.dataframe as dd


def dispatcher(_index, row):
    Row = namedtuple("row", ("npairs", "nmol", "nsample", "dist"))
    rrow = Row(*row)
    assembly = joblib.load("1A.pkl")
    try:
        result = scaffold_10x(assembly,
                              prefix="scaffold_10x", min_npairs=rrow.npairs,
                              max_dist=rrow.dist, popseq_dist=5, max_dist_orientation=5,
                              min_nsample=rrow.nsample,
                              min_nmol=rrow.nmol, unanchored=True, ncores=10, verbose=False)
        print("Test {_index}: N50".format(_index=_index), n50(result[1]["length"].values))
        joblib.dump(result, "1A_test_{}.pkl".format(_index), compress=("zlib", 6))
        return (row, result)
    except Exception:
        print("Row that failed:", row)
        raise


def grid_evaluation():

    grid = pd.DataFrame(dict(zip(("npairs", "nmol", "nsample", "dist"),
                                 list(zip(*itertools.product((2, 3),
                                                             (2, 3),
                                                             (1, 2),
                                                             range(6 * 10**4, 10**5, 10**4)))))))
    print("Starting grid evaluation")

    # pool = mp.Pool(4)
    # Row = namedtuple(, ("npairs", "nmol", "nsample", "dist"))
    # [3, 2, 2, 70000]
    # row = (3, 2, 2, 70000)
    row = (2, 2, 1, 6*10**4)

    # results = pool.starmap(dispatcher, [(0, row)])
    results = dispatcher(0, row)
                           # [(index, row) for index, row in grid.iterrows()])
    # pool.close()
    # pool.join()
    # _index, row = next(grid.iterrows())
    # result = dispatcher(assembly, _index, row)
    return results


def main():
    ne.set_num_threads(12)
    results = grid_evaluation()
    try:
        joblib.dump(results, "scaffolded.pkl", compress=("zlib", 6))
    except _pickle.PicklingError:
        print(results)
        raise


if __name__ == "__main__":
    main()
