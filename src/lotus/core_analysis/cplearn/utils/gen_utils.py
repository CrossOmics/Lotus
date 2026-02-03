import numpy as np

import hnswlib
import igraph as ig
from kneed import KneeLocator
from numba import njit


M=200
ef_construction=200
ef=200

def get_kNN(X, q=15):
    """
    Generate a k-nearest neighbors graph from the input data.
    :param X: Input data (numpy array).
    :param q: Number of nearest neighbors.
    :return: k-nearest neighbors list and distances.
    """
    n = X.shape[0]
    dim = X.shape[1]
    p = hnswlib.Index(space='l2', dim=dim)
    p.init_index(max_elements=n, ef_construction=200, M=64)
    p.add_items(X)
    p.set_ef(2*q)

    labels, dists = p.knn_query(X, k=q+1)
    knn_list = labels[:, 1:]
    knn_dists = dists[:, 1:]

    return knn_list, knn_dists


def bipartite_NN(X,Y, q=15):
    """
    Generate a k-nearest neighbors graph from the input data.
    :param X: Input data (numpy array).
    :param q: Number of nearest neighbors.
    :return: k-nearest neighbors list and distances.
    """
    n = X.shape[0]
    dim = X.shape[1]
    p = hnswlib.Index(space='l2', dim=dim)
    p.init_index(max_elements=n, ef_construction=200, M=64)
    p.add_items(X)
    p.set_ef(2*q)

    labels, dists = p.knn_query(Y, k=q+1)
    knn_list = labels[:, 1:]
    knn_dists = dists[:, 1:]

    return knn_list, knn_dists



@njit
def get_edge_list(knn_list_local,n):
    edge_list=[]
    for i_loc in range(n):
        for j_loc in knn_list_local[i_loc]:
            if i_loc != j_loc:
                edge_list.append((i_loc,j_loc))

    return edge_list



@njit
def get_induced_edge_list(edge_list,subset,n):

    hmap=-1*np.ones(n)
    for ell in subset:
        hmap[ell]=1

    new_edge_list=[]

    for u,v in edge_list:
        if hmap[u]==1 and hmap[v]==1:
            new_edge_list.append((u,v))

    return new_edge_list


def get_igraph_from_knn_list(knn_list,n):
    e_list = get_edge_list(knn_list, n)
    G = ig.Graph(directed=True)
    G.add_vertices(n)
    G.add_edges(e_list)
    G.vs["id"] = list(range(n))

    return G

from numba import types
from numba.typed import List
import numpy as np

def get_induced_knn(knn_list, knn_dist, S, n):
    """
    Filter knn_list/knn_dist so each row keeps only neighbors in S.
    Returns numba.typed.List[List[int64]], numba.typed.List[List[float64]]
    ready for @njit use.
    """
    in_subset = -1*np.ones(n, dtype=np.uint64)
    for local_id, global_id in enumerate(S):
        in_subset[global_id] = local_id


    # Ensure list-of-lists form
    if isinstance(knn_list, np.ndarray):
        knn_list = knn_list.tolist()
        knn_dist = knn_dist.tolist()

    # Top-level typed lists
    new_knn_list = List.empty_list(types.ListType(types.int64))
    new_knn_dist = List.empty_list(types.ListType(types.float64))


    for i in range(n):
        if in_subset[i]>=0:
            inner_nodes = List.empty_list(types.int64)
            inner_dists = List.empty_list(types.float64)

            nbrs=knn_list[i]
            dists=knn_dist[i]
            for j in range(len(nbrs)):

                v=int(nbrs[j])
                if in_subset[v]>=0:
                    inner_nodes.append(int(in_subset[v]))
                    inner_dists.append(float(dists[j]))


            new_knn_list.append(inner_nodes)
            new_knn_dist.append(inner_dists)

    return new_knn_list, new_knn_dist


def truncate_ng_list(ng_list, r):
    """
    Keep only the first r neighbors per node.
    Works for both list-of-lists and NumPy arrays.
    """
    # Case 1: ng_list is a NumPy array
    if isinstance(ng_list, np.ndarray):
        n, k = ng_list.shape
        if r >= k:
            return ng_list
        return ng_list[:, :r]

    # Case 2: ng_list is a Python list of lists
    out = []
    for row in ng_list:
        if len(row) > r:
            out.append(row[:r])
        else:
            out.append(row[:])
    return out



def ensure_list_of_lists(ng_list):
    """
    Converts ng_list (either np.ndarray or list of lists)
    into a numba.typed.List of lists of int64.
    """
    nb_list = List.empty_list(types.ListType(types.int64))

    # Case 1: already numpy array (n x k)
    if isinstance(ng_list, np.ndarray):
        for row in ng_list:
            inner = List.empty_list(types.int64)
            for v in row:
                inner.append(int(v))
            nb_list.append(inner)

    # Case 2: already Python list of lists
    else:
        for row in ng_list:
            inner = List.empty_list(types.int64)
            for v in row:
                inner.append(int(v))
            nb_list.append(inner)

    return nb_list


def cutoffs(values, eps=1/400, mass=0.99, win=50, rel_slope=0.01):
    v = np.asarray(values, float)
    # assume v is already sorted desc
    out = {}

    # 1) last non-zero
    nz = np.argmax(v[::-1] > 0)
    out["last_nonzero"] = len(v) - nz if v[-1] == 0 else len(v)

    # 2) absolute threshold
    i_eps = np.searchsorted(-v, -eps, side="left")
    out["value>eps"] = i_eps

    # 3) cumulative mass
    cs = np.cumsum(v)
    total = cs[-1] if cs[-1] > 0 else 1.0
    i_mass = np.searchsorted(cs, mass * total, side="left") + 1
    out["cumu_mass"] = i_mass

    # 4) knee (simple 2-line distance method)
    x = np.arange(len(values))
    y = values

    knee = KneeLocator(x, y, curve='convex', direction='decreasing')
    knee_point = knee.knee
    out["knee"] = knee_point

    # 5) flat-tail via relative slope
    vp=v[int(0.05*len(v)):]
    dv = np.abs(np.diff(vp, prepend=vp[0]))
    # rolling mean of |Îv|/v
    rel = (dv / np.maximum(vp, 1e-12))
    rel_roll = np.convolve(rel, np.ones(win)/win, mode="same")
    flat_idx = np.argmax(rel_roll < rel_slope)  # first True
    out["flat_tail"] = int(flat_idx) if rel_roll[flat_idx] < rel_slope else len(v)

    return out
