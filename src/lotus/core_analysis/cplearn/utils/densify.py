import time
import numpy as np
from numba import njit, prange

from ..utils.gen_utils import ensure_list_of_lists,truncate_ng_list
from ..utils.gen_utils import get_kNN


@njit(parallel=True)
def top_nonzero_indices_fixed(mat, f, K, subset):
    m, n = mat.shape
    out = -np.ones((m, K), dtype=np.int64)

    mat_S= mat[:, subset]

    for i in prange(m):
        # collect nonzeros
        nz_idx = []
        nz_vals = []
        for j in range(n):
            v = mat_S[i, j]
            if v != 0:
                nz_idx.append(j)
                nz_vals.append(v)

        nnz = len(nz_idx)
        if nnz == 0:
            continue

        k = f[i] if f[i] < nnz else nnz
        # sort top-k
        order = np.argsort(np.array(nz_vals))[::-1][:k]
        sel = np.array(nz_idx)[order]

        out[i, :k] = sel  # write results, rest stay -1
    return out


import numpy as np
from numba import njit, prange

@njit(parallel=True)
def topk_indices_per_row(mat_S, k):
    """
    mat: 2D array (m, n)
    k: number of top indices to extract per row
    returns: (m, k) array of indices sorted by descending values
             (ties broken by index order)
    """
    m, n = mat_S.shape

    out = np.empty((m, k), dtype=np.int64)

    for i in prange(m):
        row = mat_S[i]
        # argsort descending, then take first k
        order = np.argsort(row)[::-1][:k]
        out[i] = order
    return out



import numpy as np
from numba import njit, prange

import time


@njit(parallel=True, cache=True)
def batch_random_walk(adj_list, core_nodes, steps=5,eps=1e-4):
    """
    adj_matrix: (n, k) int array, each row has exactly k neighbors
    core_nodes: list of starting nodes
    steps: number of walk steps
    returns: (len(core_nodes), n) cumulative visitation distributions
    """
    n=len(adj_list)
    m = len(core_nodes)
    mat = np.zeros((m, n), dtype=np.float32)

    deg_inv = np.zeros(n, dtype=np.float32)
    for u in range(n):
        deg = len(adj_list[u])
        deg_inv[u] = 1.0 / deg if deg > 0 else 0.0

    max_frontier=max(1000, n//100)

    for i in prange(m):                 # parallel over seeds
        seed = core_nodes[i]

        pi = np.zeros(n, np.float32)
        visits = np.zeros(n, np.float32)
        pi[seed] = 1.0
        visits[seed] = 1.0

        new_pi = np.zeros(n, dtype=np.float32)

        # Active frontier buffers
        active = np.empty(n, np.int64)
        next_active = np.empty(n, np.int64)
        marked = np.zeros(n, np.uint8)

        active[0] = seed
        active_len = 1


        for _ in range(steps):
            new_pi.fill(0.0)
            next_len=0
            for t in range(active_len):
                u = active[t]
                p_u=pi[u]
                if p_u<=eps:
                    continue

                contrib=p_u*deg_inv[u]
                nbrs=adj_list[u]
                for j in range(len(nbrs)):
                    v = nbrs[j]
                    new_pi[v] += contrib
                    if marked[v] == 0:
                        marked[v] = 1
                        next_active[next_len] = v
                        next_len += 1


            # clear marked nodes
            for t in range(next_len):
                marked[next_active[t]] = 0




            # accumulate visits at this step
            for t in range(next_len):
                v = next_active[t]
                visits[v] += pi[v]


            # update frontier

            # ---- after constructing next_active and next_len ----
            if next_len > max_frontier:
                vals = np.empty(next_len, np.float32)
                for t in range(next_len):
                    vals[t] = new_pi[next_active[t]]

                # find threshold for top `max_frontier` nodes
                k_drop = next_len - max_frontier
                thr = np.partition(vals, k_drop)[k_drop]

                # retain only high-mass nodes
                new_len = 0
                for t in range(next_len):
                    v = next_active[t]
                    if new_pi[v] >= thr:
                        next_active[new_len] = v
                        new_len += 1
                next_len = new_len

            #Swap current distributions
            pi, new_pi = new_pi, pi

            #Now swap active and next_active
            active, next_active = next_active, active
            active_len = next_len


            # for t in range(next_len):
            #     active[t] = next_active[t]
            # active_len = next_len

            if active_len == 0:
                break


        mat[i] = visits

    return mat


def degree_within_subset_fast(adj_list, subset):
    n = len(adj_list)
    mask = np.zeros(n, dtype=bool)
    mask[subset] = True

    deg_full = np.zeros(len(subset), dtype=int)
    deg_sub = np.zeros(len(subset), dtype=int)

    for i, u in enumerate(subset):
        nbrs = np.array(adj_list[u]).astype(np.int32)
        deg_full[i] = len(nbrs)
        deg_sub[i] = np.count_nonzero(mask[nbrs])
    return deg_full, deg_sub



def densify_rw(adj_list,subset_nodes):

    t0 = time.time()

    n=len(adj_list)
    n_core=len(subset_nodes)

    deglist_full,deglist_sub = degree_within_subset_fast(adj_list, subset_nodes)

    hmap=-1*np.ones(n, dtype=np.int64)
    for i,u in enumerate(subset_nodes):
        hmap[u]=i


    set_to_expand=[]
    vacancy=[]

    total_vacancy=0
    for i in range(n_core):
        if deglist_sub[i]<0.75*deglist_full[i]:


            set_to_expand.append(int(subset_nodes[i]))
            vacancy.append(int(deglist_full[i]-deglist_sub[i]))
            total_vacancy += int(deglist_full[i]-deglist_sub[i])



    set_to_expand=np.array(set_to_expand).astype(int)

    print("Number of nodes to expand and total number edges that could be added:",len(set_to_expand)," ",total_vacancy,flush=True)


    if len(set_to_expand) == 0 or max(vacancy) == 0:
        return adj_list

    ste_map=-1*np.ones(n,dtype=int)
    for i,u in enumerate(set_to_expand):
        ste_map[u]=i


    t1=time.time()

    adj_list_numba=ensure_list_of_lists(adj_list)

    mat=batch_random_walk(adj_list_numba, set_to_expand)


    t2=time.time()

    #Candidate_nbrs shape is (len(set_to_expand), max(vacancy))
    mat_S= mat[:, np.array(subset_nodes).astype(int)]
    candidate_nbrs=topk_indices_per_row(mat_S, max(vacancy))

    t3=time.time()

    #Create densified adjacency list. Just send whole graph.
    adj_list_dense=[list() for _ in range(n)]

    for i in range(n):

        if ste_map[i]==-1:
            adj_list_dense[i]=adj_list[i]

        else:
            nbrs=[]
            for v in adj_list[i]:
                    nbrs.append(v)

            c_temp=0
            for v in candidate_nbrs[ste_map[i]]:
                if mat[ste_map[i],v]==0:
                    break
                if hmap[v]!=-1:
                    nbrs.append(v)
                    c_temp+=1

                if c_temp==vacancy[ste_map[i]]:
                    break

            adj_list_dense[i] = nbrs

        t4=time.time()


    print(f"Time for step 1: {t1-t0:.3f}")
    print(f"Time for step 2: {t2-t1:.3f}")
    print(f"Time for step 3: {t3-t2:.3f}")
    print(f"Time for step 4: {t4-t3:.3f}")

    return adj_list_dense


def densify_knn(X,adj_list,subset_nodes):
    n=len(adj_list)
    r=max(len(row) for row in adj_list)
    print("Original average degree:", np.mean([len(row) for row in adj_list]), flush=True)


    knn_list_sub, knn_dist_sub = get_kNN(X[subset_nodes], r)

    adj_list_dense=[list() for _ in range(n)]

    hmap=-1*np.ones(n, dtype=np.int64)
    for i,u in enumerate(subset_nodes):
        hmap[u]=i

    for i in range(n):

        if hmap[i]==-1:
            adj_list_dense[i]=adj_list[i]

        else:
            nbrs=[]
            for v in adj_list[i]:
                    nbrs.append(v)

            for v in knn_list_sub[hmap[i]]:
                if hmap[v]!=-1 and v not in nbrs:
                    nbrs.append(v)

            adj_list_dense[i] = nbrs





    return adj_list_dense



