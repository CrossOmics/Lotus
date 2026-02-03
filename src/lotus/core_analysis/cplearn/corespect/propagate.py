
from.diffusion import label_diffusion
import numpy as np
import scipy.sparse as sp

def adjlist_to_csr(adj_list):
    n = len(adj_list)
    rows, cols = [], []
    for u, nbrs in enumerate(adj_list):
        rows.extend([u] * len(nbrs))
        cols.extend(nbrs)
    data = np.ones(len(rows), dtype=np.float32)
    return sp.csr_matrix((data, (rows, cols)), shape=(n, n))

import numpy as np

def divide_into_layers(F, P, x=5, method='argmax_conf'):
    """
    Divide non-core points into x layers per class based on confidence.

    Parameters
    ----------
    F : np.ndarray (n_samples, n_classes)
        Probability matrix from label propagation or classifier.
    P : array-like
        Indices of non-core points.
    x : int
        Number of layers.
    method : {'argmax_conf', 'margin', 'entropy'}
        Scoring method for ranking points:
        - 'argmax_conf': highest probability (standard confidence)
        - 'margin': difference between top two probabilities
        - 'entropy': lower entropy = higher confidence
    """
    Fp = F[P]
    n_classes = F.shape[1]

    # Compute scores based on method
    if method == 'argmax_conf':
        scores = Fp.max(axis=1)
    elif method == 'margin':
        part_sorted = np.partition(Fp, -2, axis=1)
        scores = part_sorted[:, -1] - part_sorted[:, -2]
    elif method == 'entropy':
        eps = 1e-12
        scores = -np.sum(Fp * np.log(Fp + eps), axis=1)
        scores = -scores  # higher = more confident
    else:
        raise ValueError("Unknown method")

    labels = Fp.argmax(axis=1)
    layers = [[] for _ in range(x)]

    # For each class, split its points into x equal-sized bins
    for c in range(n_classes):
        idx_c = np.where(labels == c)[0]
        if len(idx_c) == 0:
            continue
        order = np.argsort(-scores[idx_c])  # descending confidence
        split_points = np.array_split(idx_c[order], x)
        for i, subset in enumerate(split_points):
            layers[i].extend(P[subset])  # map back to original indices

    # Convert to np arrays
    layers = [np.array(layer, dtype=int) for layer in layers]
    return layers



def propagate_from_core(adj_list,layers_,labels_,connectivity,normalized,mode,alpha,max_iter,**kwargs):


    n=len(adj_list)
    curr_nodes=[]

    new_layers_=[]

    #Changing rem_nodes to be everything but the core nodes temporarily.
    #ToDo: Add an option.
    for layer in layers_[0:1]:
        curr_nodes.extend(layer)
        new_layers_.append(layer)


    rem_nodes=np.array(list(set(range(n)) - set(curr_nodes))).astype(int)

    if connectivity=='out':
        A_csr = adjlist_to_csr(adj_list)

    elif connectivity=='in':
        A_csr = adjlist_to_csr(adj_list).T

    elif connectivity=='sym':
        A_csr = adjlist_to_csr(adj_list)
        A_csr = 0.5*(A_csr + A_csr.T)

    else:
        raise ValueError(f"Invalid direction: {connectivity}")


    if normalized:
        norm_mode='normalized'
    else:
        norm_mode='random_walk'

    if mode == 'default':
        alpha=1.0



    labels_, F, F_history = label_diffusion(
        A_csr,
        labels_,
        norm_mode=norm_mode,
        alpha=alpha,
        max_iter=max_iter,
        tol=1e-5,
        dtype=np.float32,
        num_threads=8,
        save_every=5,
    )

    additional_layers=divide_into_layers(F,rem_nodes,x=5,method='argmax_conf')

    for layer in additional_layers:
        new_layers_.append(layer)




    return labels_,new_layers_,F