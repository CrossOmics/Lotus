import numpy as np

from joblib import Parallel, delayed

from ..utils.gen_utils import ensure_list_of_lists
from ..utils.cluster_utils import cluster_subset

from numba import njit, prange


@njit(parallel=True)
def build_edges_and_out_degree_numba(adj_list):
    n = len(adj_list)

    # Count total edges
    # FIXME: precompute the offset for each vertex
    start_indices = np.zeros(n, dtype=np.int64)
    total_edges = 0
    for u in range(n):
        start_indices[u] = total_edges
        total_edges += len(adj_list[u])

    edge_list = np.empty((total_edges, 2), dtype=np.int64)
    out_degree = np.empty(n, dtype=np.int64)

    # FIXME: refactored to avoid race condition on the global variable `offset`
    # Problem: using a global `offset` in a parallel loop results in race conditions, leading to IndexError
    # Solution: use a precomputed `start_indices` array to store an independent offset for each vertex
    # This is equivalent to the original version.

    for u in prange(n):
        u_offset = start_indices[u]
        nbrs = adj_list[u]
        out_degree[u] = len(nbrs)
        for i in range(len(nbrs)):
            edge_list[u_offset + i, 0] = u
            edge_list[u_offset + i, 1] = nbrs[i]

    return edge_list, out_degree


def majority_single_iteration(edge_list, deg_list, n, cnum, init_labels, rem_nodes_mask, eps=0):
    votes = np.zeros((n, cnum))
    added_nodes = []
    for i in range(len(edge_list)):
        u, v = edge_list[i]
        # for u, v in edge_list:
        if v >= len(init_labels):
            print(f"\ninvalid index! v={v}, n={len(init_labels)}, edge index={i}")
            continue
        c_idx = int(init_labels[v])
        if rem_nodes_mask[u] == 1 and c_idx != -1:
            votes[u, c_idx] += 1

    for u in range(n):

        if rem_nodes_mask[u] == 1:
            if max(votes[u]) > (0.5 + eps) * deg_list[u]:
                init_labels[u] = np.argmax(votes[u])
                rem_nodes_mask[u] = 0
                added_nodes.append(u)

    return init_labels, rem_nodes_mask, added_nodes


def recursive_majority(adj_list_n, core_nodes, rem_nodes, resolution=1):
    adj_list = ensure_list_of_lists(adj_list_n)

    n = len(adj_list)
    tolerance = int(0.005 * n)
    max_iter = 500  # SafeGuard to avoid infinite loops

    init_labels = cluster_subset(adj_list, core_nodes, resolution=resolution)

    cnum = len(set(init_labels)) - (1 if -1 in init_labels else 0)
    edge_list, deg_list = build_edges_and_out_degree_numba(adj_list)

    rem_nodes_mask = np.zeros(n).astype(int)
    rem_nodes_mask[rem_nodes] = 1

    layer = []

    round_num = 0
    while True:
        n0 = np.sum(init_labels != -1)
        proxy_labels, rem_nodes_mask, added_nodes = majority_single_iteration(edge_list, deg_list, n, cnum, init_labels,
                                                                              rem_nodes_mask)

        for ell in added_nodes:
            layer.append(ell)

        n1 = np.sum(proxy_labels != -1)

        if n1 - n0 < tolerance:
            return layer

        if n0 > 0.95 * n:
            return layer

        if round_num > max_iter:
            return layer

        round_num += 1


def find_intersection_layers(adj_list, core_nodes, rem_nodes, ng_num, resolution):
    n = len(adj_list)

    if rem_nodes is None:
        rem_nodes = np.array([i for i in range(n) if i not in core_nodes])

    np.random.shuffle(rem_nodes)  # Should not effect reproducibility

    layer_list = Parallel(n_jobs=10)(
        delayed(recursive_majority)(adj_list, core_nodes, rem_nodes, resolution) for _ in range(10))

    layer_candidate = set(layer_list[0])
    for i in range(1, 10):
        layer_candidate = layer_candidate.intersection(set(layer_list[i]))

    return layer_candidate


def choose_stopping_res(adj_list, core_nodes, thr=0.95, starting_resolution=1):
    n0 = len(core_nodes)

    n1 = n0
    t = starting_resolution
    while n1 > thr * n0 and t < 5:
        layer_candidate = find_intersection_layers(adj_list, core_nodes, ng_num=20,
                                                   rem_nodes=np.array(core_nodes).astype(int),
                                                   resolution=starting_resolution)
        n1 = len(layer_candidate)
        print(t, n1, n0)
        t += 0.25

    return t - 0.5
