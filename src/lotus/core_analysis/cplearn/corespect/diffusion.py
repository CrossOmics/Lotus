import os
import numpy as np
import scipy.sparse as sp

def label_diffusion(
    A_csr: sp.csr_matrix,
    y: np.ndarray,
    *,
    norm_mode: str = "normalized",   # "random_walk" -> D^{-1}A, "normalized" -> D^{-1/2}AD^{-1/2}
    alpha: float = 0.9,
    max_iter: int = 50,
    tol: float = 1e-5,
    dtype=np.float32,
    num_threads: int = 8,
    save_every: int = 5,              # <-- snapshot cadence
):
    """
    Diffusion-based label propagation with absorbing cores.
    Saves F every `save_every` iterations, plus t=0 and the final state.

    Returns
    -------
    hard_labels : (n,) int
    F           : (n, K) float
    F_history   : (T_saves, n, K) float  (snapshots)
    """
    # threading hints
    os.environ["OPENBLAS_NUM_THREADS"] = str(num_threads)
    os.environ["OMP_NUM_THREADS"] = str(num_threads)

    if not sp.isspmatrix_csr(A_csr):
        A_csr = A_csr.tocsr(copy=False)

    n = A_csr.shape[0]
    y = np.asarray(y, dtype=np.int64)
    if y.shape[0] != n:
        raise ValueError("y must have length equal to A_csr.shape[0]")

    core_mask = (y >= 0)
    if not np.any(core_mask):
        return y.copy(), np.zeros((n, 0), dtype=dtype), np.zeros((0, n, 0), dtype=dtype)

    classes = np.unique(y[core_mask])
    K = len(classes)
    lab2idx = {c: i for i, c in enumerate(classes)}

    # ---- build diffusion operator P (O(nnz))
    A_csr = A_csr.astype(dtype, copy=False)
    d = np.asarray(A_csr.sum(axis=1)).ravel().astype(dtype)
    d_safe = np.maximum(d, np.finfo(dtype).tiny)

    P = A_csr.copy()
    indptr, indices, data = P.indptr, P.indices, P.data
    if norm_mode == "random_walk":
        invdeg = 1.0 / d_safe
        data *= np.repeat(invdeg, np.diff(indptr))
    elif norm_mode == "normalized":
        inv_sqrt = 1.0 / np.sqrt(d_safe)
        for i in range(n):
            s, e = indptr[i], indptr[i+1]
            if e > s:
                data[s:e] *= inv_sqrt[i] * inv_sqrt[indices[s:e]]
    else:
        raise ValueError("norm_mode must be 'random_walk' or 'normalized'")
    P.eliminate_zeros()

    # ---- initialize F
    F = np.zeros((n, K), dtype=dtype)
    for i, lab in enumerate(y):
        if lab >= 0:
            F[i, lab2idx[lab]] = 1.0
    F_core = F[core_mask].copy()

    # ---- history setup
    F_history = []
    last_saved_t = -1
    if save_every and save_every > 0:
        F_history.append(F.copy())  # t=0
        last_saved_t = 0

    # ---- diffusion loop
    for t in range(1, max_iter + 1):
        F_new = (1 - alpha) * F + alpha * (P @ F)
        F_new[core_mask] = F_core

        # snapshot every `save_every` iters
        if save_every and save_every > 0 and (t % save_every == 0):
            F_history.append(F_new.copy())
            last_saved_t = t

        if np.max(np.abs(F_new - F)) < tol:
            F = F_new
            break
        F = F_new
    else:
        # exhausted all iterations without tol break
        pass

    # ensure final snapshot is included
    if save_every and save_every > 0 and last_saved_t != t:
        F_history.append(F.copy())

    # ---- hard labels
    row_max = F.max(axis=1)
    argmax = F.argmax(axis=1)
    hard_labels = np.full(n, -1, dtype=int)
    nonzero = row_max > 0
    hard_labels[nonzero] = classes[argmax[nonzero]]

    # stack history -> (T_saves, n, K)
    F_history = np.stack(F_history, axis=0) if F_history else np.zeros((0, n, K), dtype=dtype)
    return hard_labels, F, F_history
