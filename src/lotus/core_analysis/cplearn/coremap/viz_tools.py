from sympy.physics.control import Parallel

from .find_anchors import find_anchors
from .coremap_graph_utils import varying_heat_kernel
from ..utils.gen_utils import bipartite_NN

import numpy as np
from numba import njit
from numba import prange

def _make_epochs_per_sample(weights, n_epochs):
    w = np.asarray(weights, dtype=np.float64)
    if w.size == 0:
        return w
    result = np.full(w.shape[0], -1.0, dtype=np.float64)  # sentinel: never sample
    wmax = float(w.max())
    if wmax <= 0.0:
        return result
    n_samples = n_epochs * (w / wmax)
    mask = n_samples > 0.0
    result[mask] = float(n_epochs) / n_samples[mask]
    return result


import numpy as np


def get_adj_dists(adj_list, X, metric="euclidean"):
    """
    Sorts each adjacency list by increasing distance between node i and its neighbors.

    Parameters
    ----------
    adj_list : list[list[int]]
        Adjacency list (neighbors for each node)
    X : np.ndarray, shape (n, d)
        Node feature matrix
    metric : str, optional
        Distance metric ('euclidean' or 'cosine')

    Returns
    -------
    sorted_adj : list[list[int]]
        Sorted adjacency list by distance
    sorted_dists : list[list[float]]
        Sorted list of corresponding distances
    """
    n = len(adj_list)
    sorted_adj = []
    sorted_dists = []

    for i in range(n):
        neighbors = np.array(adj_list[i]).astype(int)
        if neighbors.size == 0:
            sorted_adj.append([])
            sorted_dists.append([])
            continue

        # Compute distances to all neighbors
        if metric == "euclidean":
            dists = np.linalg.norm(X[neighbors] - X[i], axis=1)
        elif metric == "cosine":
            xi = X[i]
            sims = X[neighbors] @ xi / (
                    np.linalg.norm(X[neighbors], axis=1) * np.linalg.norm(xi) + 1e-12
            )
            dists = 1 - sims
        else:
            raise ValueError("metric must be 'euclidean' or 'cosine'")

        # Sort neighbors by distance
        order = np.argsort(dists)
        sorted_adj.append(neighbors[order].tolist())
        sorted_dists.append(dists[order].tolist())

    return sorted_adj, sorted_dists


def _create_anchored_edge_list(corespect_obj,anchor_list,anchor_distances,X_umap,anchor_finding_mode,final_prob_vec,anchor_reach=None):


    core_layer=corespect_obj.layers_[0]

    X_p=corespect_obj.X
    knn_list=corespect_obj.adj_list
    knn_list,knn_dists=get_adj_dists(knn_list,X_p,metric='euclidean')


    P_vec = varying_heat_kernel(knn_dists)

    #Get the normal umap-weighted graph
    edges_from=[]
    edges_to=[]
    weights=[]
    for u in range(len(knn_list)):
        for j,v in enumerate(knn_list[u]):
            edges_from.append(u)
            edges_to.append(v)
            weights.append(P_vec[u][j])


    #Get the anchored graph

    anchor_coordinates=[]
    for u in range(len(anchor_list)):
        anchor_coordinates.append(np.mean(X_umap[anchor_list[u]],axis=0))

    if anchor_reach is None:
        anchor_reach=0.1*len(core_layer)

    if anchor_finding_mode == 'default':
        P_vec_anchor=varying_heat_kernel(anchor_distances,th=anchor_reach)

        anchored_edges_from=[]
        anchor_to=[]
        anchored_weights=[]
        for u in range(len(anchor_list)):
            for j,v in enumerate(anchor_list[u]):
                anchored_edges_from.append(v)
                anchor_to.append(anchor_coordinates[u])
                anchored_weights.append(P_vec_anchor[u][j])


    elif anchor_finding_mode == 'smoothed':
        fx=0
        P_vec_anchor = varying_heat_kernel(anchor_distances, th=anchor_reach)

        anchored_edges_from = []
        anchor_to = []
        anchored_weights = []
        for u in range(len(anchor_list)):
            for j, v in enumerate(anchor_list[u]):
                anchored_edges_from.append(v)
                anchor_to.append(anchor_coordinates[u])

                if u in final_prob_vec[v].keys():
                    anchored_weights.append(P_vec_anchor[u][j]*final_prob_vec[v][u]) #Normalizing with GMM_probability.
                    fx+=1

                else:
                    anchored_weights.append(P_vec_anchor[u][j])


        print("Total weights distributed=",fx)

    else:
        raise KeyError('Unknown anchor_finding_mode')



    return edges_from,edges_to,weights,anchored_edges_from,anchor_to,anchored_weights

@njit
def clip(val):
    if val > 4.0:
        return 4.0
    elif val < -4.0:
        return -4.0
    else:
        return val


@njit
def clip_th(val,th=4.0):
    if val > th:
        return th
    elif val < -th:
        return -th
    else:
        return val


@njit(parallel=True)
def _optimize_layout_euclidean_single_epoch(
    coordinates,            # (n, dim)
    weights,                # (|E|,)
    edge_from,              # (|E|,)
    edge_to,                # (|E|,)
    epochs_per_sample,      # (|E|,)
    epoch_of_next_sample,   # (|E|,)
    epochs_per_negative_sample,        # (|E|,)
    epoch_of_next_negative_sample,     # (|E|,)
    n,                      # current epoch index
    alpha,                  # learning rate
    a, b,
    repulsion_ratio
):
    n_vertices = coordinates.shape[0]
    dim = coordinates.shape[1]

    #Added prange.
    for i in prange(epochs_per_sample.shape[0]):

        # ---- Positive (attractive) sample ----
        if epoch_of_next_sample[i] <= n:
            j = edge_from[i]
            k = edge_to[i]

            current = coordinates[j]
            other   = coordinates[k]

            # squared distance
            dist2 = 0.0
            for d in range(dim):
                diff = current[d] - other[d]
                dist2 += diff * diff

            if dist2 > 0.0:
                grad_coeff = -2.0 * a * b * (dist2 ** (b - 1.0))
                grad_coeff /= (a * dist2 + 1.0)
            else:
                grad_coeff = 0.0


            # Update BOTH endpoints for attraction
            for d in range(dim):
                grad_d = clip(grad_coeff * (current[d] - other[d]))
                delta  = grad_d * alpha
                current[d] += delta
                other[d]   -= delta  # symmetric update

            epoch_of_next_sample[i] += epochs_per_sample[i]

            # ---- Negative (repulsive) samples ----
            # Count how many negative updates are due at epoch n
            n_neg = 0
            if epochs_per_negative_sample[i] < 1e308:  # not inf
                if epoch_of_next_negative_sample[i] <= n:
                    # include current epoch; floor division equivalent
                    n_neg = int((n - epoch_of_next_negative_sample[i]) /
                                epochs_per_negative_sample[i]) + 1
                    if n_neg < 0:
                        n_neg = 0  # safety

            for _ in range(n_neg):
                # simple RNG; reproducibility not guaranteed under njit
                kk = np.random.randint(n_vertices)
                if kk == j:
                    continue
                other = coordinates[kk]

                dist2 = 0.0
                for d in range(dim):
                    diff = current[d] - other[d]
                    dist2 += diff * diff

                if dist2 > 0.0:
                    grad_coeff = (2.0 * b) / ((0.001 + dist2) * (a * (dist2 ** b) + 1.0))
                else:
                    grad_coeff = 0.0

                if grad_coeff > 0.0:
                    for d in range(dim):
                        grad_d = clip(grad_coeff * (current[d] - other[d]))
                        current[d] += grad_d * alpha

            epoch_of_next_negative_sample[i] += n_neg * epochs_per_negative_sample[i]




@njit(parallel=True)
def _stabilize_layout_with_anchor_single_epoch(
    coordinates,            # (n, dim),            # (n, dim)
    weights,                # (|E|,)
    edge_from,              # (|E|,)
    coordinates_to,                # (|E|,)
    epochs_per_sample,      # (|E|,)
    epoch_of_next_sample,   # (|E|,)
    n,                      # current epoch index
    alpha,                  # learning rate
    a, b
):
    dim = coordinates.shape[1]


    #Added prange.
    for i in prange(epochs_per_sample.shape[0]):

        # ---- Positive (attractive) sample ----
        if epoch_of_next_sample[i] <= n:
            j = edge_from[i]
            current = coordinates[j]
            other=coordinates_to[i]

            #print("check indices",current,other)

            # squared distance
            dist2 = 0.0
            for d in range(dim):
                diff = current[d] - other[d]
                dist2 += diff * diff

            if dist2 > 0.0:
                grad_coeff = -2.0 * a * b * (dist2 ** (b - 1.0))
                grad_coeff /= (a * dist2 + 1.0)
            else:
                grad_coeff = 0.0


            # Update BOTH endpoints for attraction
            for d in range(dim):
                grad_d = clip(grad_coeff * (current[d] - other[d]))
                delta  = grad_d * alpha
                current[d] += delta
                #other[d]   -= delta

            epoch_of_next_sample[i] += epochs_per_sample[i]



def anchored_map_single_layer(edge_from_global,edge_to_global,weights_global,anchored_edge_from_global,anchor_to_global,anchored_weights_global,init_coordinate,layer,total_size):


    hmap=-1*np.ones(total_size)
    for i,node in enumerate(layer):
        hmap[node]=i


    edge_from=[]
    edge_to=[]
    weights=[]

    for ell in range(len(edge_from_global)):

        if hmap[edge_from_global[ell]]>=0 and hmap[edge_to_global[ell]]>=0:

            edge_from.append(hmap[edge_from_global[ell]])
            edge_to.append(hmap[edge_to_global[ell]])
            weights.append(weights_global[ell])


    edge_from = np.asarray(edge_from, dtype=np.int64)
    edge_to = np.asarray(edge_to, dtype=np.int64)
    weights = np.asarray(weights, dtype=np.float64)

    anchored_edge_from=[]
    anchor_to=[]
    anchored_weights=[]

    # Setting umap parameters
    total_epochs = 400
    repulsion_ratio = 5.0
    initial_alpha = 1.0
    a = 1.579
    b = 0.895

    epochs_per_sample = _make_epochs_per_sample(weights, total_epochs)

    # schedule: never-sampled edges => +inf, not negative numbers
    epoch_of_next_sample = np.where(epochs_per_sample > 0.0,
                                    epochs_per_sample.copy(),
                                    np.inf)

    epochs_per_negative_sample = np.where(epochs_per_sample > 0.0,
                                          epochs_per_sample / repulsion_ratio,
                                          np.inf)
    epoch_of_next_negative_sample = epochs_per_negative_sample.copy()


    for ell in range(len(anchored_edge_from_global)):
        if hmap[anchored_edge_from_global[ell]]>=0:
            anchored_edge_from.append(hmap[anchored_edge_from_global[ell]])
            anchor_to.append(anchor_to_global[ell])
            anchored_weights.append(anchored_weights_global[ell])


    anchor_edge_from = np.asarray(anchored_edge_from, dtype=np.int64)
    anchor_to = np.asarray(anchor_to, dtype=np.float64)
    anchor_weights = np.asarray(anchored_weights, dtype=np.float64)

    #print(np.shape(anchor_edge_from), np.shape(anchor_to), np.shape(anchor_weights))

    anchor_epochs_per_sample = _make_epochs_per_sample(anchor_weights, total_epochs)

    # schedule: never-sampled edges => +inf, not negative numbers
    anchor_epoch_of_next_sample = np.where(anchor_epochs_per_sample > 0.0,
                                           anchor_epochs_per_sample.copy(),
                                           np.inf)

    alpha = float(initial_alpha)

    coordinates=init_coordinate

    for n in range(total_epochs):
        _optimize_layout_euclidean_single_epoch(
            coordinates,
            weights, edge_from, edge_to,
            epochs_per_sample,
            epoch_of_next_sample,
            epochs_per_negative_sample,
            epoch_of_next_negative_sample,
            n, alpha, a, b, repulsion_ratio
        )
        alpha = initial_alpha * (1.0 - (n + 1.0) / float(total_epochs))
        if alpha < 0.0:
            alpha = 0.0

        _stabilize_layout_with_anchor_single_epoch(
            coordinates,
            anchor_weights, anchor_edge_from, anchor_to,
            anchor_epochs_per_sample,
            anchor_epoch_of_next_sample,
            n, alpha, a, b
        )

    return coordinates


#The anchors should be centroid?
#Do I move the initial position of the anchors around?


def anchored_map(coremap):

    anchor_list,anchor_distances,final_prob_vec=find_anchors(coremap.core_obj,anchor_finding_mode=coremap.anchor_finding_mode)

    edge_from,edge_to,weights,anchored_edge_from,anchor_to,anchored_weights=_create_anchored_edge_list(coremap.core_obj,anchor_list,anchor_distances,coremap.X_umap,coremap.anchor_finding_mode,final_prob_vec,coremap.anchor_reach)



    Layers=coremap.layers_

    prev_layer=None
    prev_coordinate=None

    label_dict={}


    curr_layer=[]

    for round,layer in enumerate(Layers):

        if prev_layer is None:
            coordinate=coremap.X_umap[layer]
            curr_layer.extend(layer)

        else:
            curr_layer.extend(layer)
            coordinate = coremap.X_umap[curr_layer]

            prev_idx=np.array([i for i in range(len(prev_layer))]).astype(int)
            coordinate[prev_idx]=prev_coordinate


        print(f"Shape of embedding after round {round} is {coordinate.shape}")

        coordinate=anchored_map_single_layer(edge_from,edge_to,weights,anchored_edge_from,anchor_to,anchored_weights,coordinate,curr_layer,coremap.core_obj.X.shape[0])
        prev_coordinate=coordinate.copy()
        prev_layer=curr_layer.copy()

        label_dict[round]=coordinate



    return label_dict



