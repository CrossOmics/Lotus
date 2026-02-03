import numpy as np
from numba import njit
from numba import prange

@njit
def find_sigma(distances,rhos,i, target_i):
    lo, hi = 1e-8, 100000.0  # Search range for sigma
    for _ in range(64):  # Binary search
        sigma = (lo + hi) / 2.0
        weights = np.exp(-(distances[i] - rhos[i]) / sigma)
        if np.sum(weights) > target_i:
            hi = sigma
        else:
            lo = sigma
    return (lo + hi) / 2.0

def varying_heat_kernel(distances,return_mode='default',th=None):
    n_points = len(distances)
    sigmas = np.zeros(n_points)
    target = np.zeros(n_points)
    rhos = np.zeros(n_points)

    P_vec = [[] for _ in range(n_points)]

    # Step 1: Compute rho_i (local connectivity)
    for i in range(n_points):

        rhos[i] = distances[i][0] if len(distances[i]) > 0 else 0  # Minimum nonzero distance

        if rhos[i] != min(distances[i]):
            raise KeyError(f"0-th index is not closest, {rhos[i]:.3f} {min(distances[i]):.3f}")

    # Step 2: Solve for sigma_i using binary search


    for i in prange(n_points):
        if len(distances[i]) > 0:
            c1 = len(distances[i])

            if c1 == 1:
                target[i] = 1

            elif c1 <= 5:
                target[i] = 1 + c1 / 2

            else:
                if th is None:
                    target[i] = np.log2(c1)
                else:
                    target[i] = th

            #in_place
            lo, hi = 1e-8, 100000.0  # Search range for sigma
            for _ in range(64):  # Binary search
                sigma = (lo + hi) / 2.0
                weights = np.exp(-(distances[i] - rhos[i]) / sigma)
                if np.sum(weights) > target[i]:
                    hi = sigma
                else:
                    lo = sigma

            sigmas[i] = (lo + hi) / 2.0

    # Step 3: Compute edge weights
    for i in prange(n_points):
        for d in distances[i]:
            weight_vec = np.exp(-(d - rhos[i]) / sigmas[i])

            P_vec[i].append(weight_vec)


    if return_mode=='anchor':
        return P_vec, sigmas, rhos, target

    return P_vec


import networkx as nx






