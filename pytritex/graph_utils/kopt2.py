import pandas as pd
import numpy as np


# import numpy as np
#
# # Calculate the euclidian distance in n-space of the route r traversing cities c, ending at the path start.
# path_distance = lambda r,c: np.sum([np.linalg.norm(c[r[p]]-c[r[p-1]]) for p in range(len(r))])
# # Reverse the order of all elements from element i to element k in array r.
# two_opt_swap = lambda r,i,k: np.concatenate((r[0:i],r[k:-len(r)+i-1:-1],r[k+1:len(r)]))
#
# def two_opt(cities,improvement_threshold): # 2-opt Algorithm adapted from https://en.wikipedia.org/wiki/2-opt
#     route = np.arange(cities.shape[0]) # Make an array of row numbers corresponding to cities.
#     improvement_factor = 1 # Initialize the improvement factor.
#     best_distance = path_distance(route,cities) # Calculate the distance of the initial path.
#     while improvement_factor > improvement_threshold: # If the route is still improving, keep going!
#         distance_to_beat = best_distance # Record the distance at the beginning of the loop.
#         for swap_first in range(1,len(route)-2): # From each city except the first and last,
#             for swap_last in range(swap_first+1,len(route)): # to each of the cities following,
#                 new_route = two_opt_swap(route,swap_first,swap_last) # try reversing the order of these cities
#                 new_distance = path_distance(new_route,cities) # and check the total distance with this modification.
#                 if new_distance < best_distance: # If the path distance is an improvement,
#                     route = new_route # make this the accepted best route
#                     best_distance = new_distance # and update the distance corresponding to this route.
#         improvement_factor = 1 - best_distance/distance_to_beat # Calculate how much the route has improved.
#     return route # When the route is no longer improving substantially, stop searching and return the route.


def checker(df, edge_list):
    # df has the columns "cluster", "bin"
    m = pd.DataFrame({"cluster1": df["cluster"][:-1].values, "cluster2": df["cluster"][1:].values})
    # Edge list is indexed on "cluster1", "cluster2"
    m = edge_list.merge(m, left_index=True, right_on=["cluster1", "cluster2"], how="right")[
        ["cluster1", "cluster2", "weight"]]
    m.columns = ["cluster1", "cluster2", "weight12"]
    __left = df.copy()  # "cluster", "bin"
    __left.columns = ["cluster1", "bin1"]
    __left = __left.set_index("cluster1")
    m = __left.merge(m.set_index("cluster1").sort_index(axis=1),
                     right_index=True, left_index=True).reset_index(drop=False)
    __left = df.copy()
    __left.columns = ["cluster2", "bin2"]
    __left = __left.set_index("cluster2")
    m = __left.merge(m.set_index("cluster2").sort_index(axis=1),
                     how="right", left_index=True, right_index=True).reset_index(drop=False)
    n = m.copy().rename(columns={"cluster1": "cluster3",
                                 "cluster2": "cluster4",
                                 "bin1": "bin3",
                                 "bin2": "bin4",
                                 "weight12": "weight34"}, errors="raise")
    # Cartesian product using a dummy key.
    mn = n.assign(dummy=1).set_index("dummy").sort_index(axis=1).merge(
        m.assign(dummy=1).set_index("dummy").sort_index(axis=1),
        left_index=True, right_index=True, how="outer").query("bin1 < bin3").reset_index(drop=True)
    # Cluster1, cluster2 are in the index
    o = edge_list.loc[:, "weight"].reset_index(drop=False)
    key = ["cluster1", "cluster3"]
    __left = o.copy()
    __left.columns = key + ["weight13"]
    mn = __left.set_index(key).merge(mn.set_index(key),
                                     left_index=True,
                                     right_index=True).reset_index(drop=False)
    key = ["cluster2", "cluster4"]
    __left = o.copy()
    __left.columns = key + ["weight24"]
    mn = __left.set_index(key).merge(mn.set_index(key), left_index=True, right_index=True).reset_index(drop=False)
    mn.eval("old = weight12 + weight34", inplace=True)
    mn.eval("new = weight13 + weight24", inplace=True)
    mn.eval("diff = old - new", inplace=True)
    mn = mn.sort_values("diff", ascending=False)
    return mn


# Traveling salesman heuristic for Hi-C mapping construction
def kopt2(df: pd.DataFrame, edge_list: pd.DataFrame):
    # df has the columns "cluster", "bin"
    edge_list = edge_list.set_index(["cluster1", "cluster2"])
    mn = checker(df, edge_list)
    while mn.loc[mn["diff"] > 0].shape[0] > 0:
        x = mn.head(1)
        bin1 = df.query("cluster == @x.cluster1")["bin"]
        bin2 = df.query("cluster == @x.cluster2")["bin"]
        bin3 = df.query("cluster == @x.cluster3")["bin"]
        bin4 = df.query("cluster == @x.cluster4")["bin"]
        # bin1 = df.loc[df["cluster"] == x["cluster1"], "bin"]
        # bin2 = df.loc[df["cluster"] == x["cluster2"], "bin"]
        # bin3 = df.loc[df["cluster"] == x["cluster3"], "bin"]
        # bin4 = df.loc[df["cluster"] == x["cluster4"], "bin"]
        df = pd.concat([df.iloc[:bin1], df.iloc[bin3:bin2], df.iloc[bin4:]]).reset_index(drop=True)
        df = df.assign(bin=np.arange(df.shape[0], dtype=np.int))
        mn = checker(df, edge_list)
    return df
