
from ..utils.densify import densify_knn
from ..utils.stable_core_utils import choose_stopping_res

import leidenalg
import numpy as np

import time

import igraph as ig


from sklearn.mixture import GaussianMixture

from tqdm import tqdm


def should_continue(prev_bic, curr_bic, delta_bic_min=10.0):
    # continue if improvement is at least -delta_bic_min (e.g., -10 or lower)

    delta_bic_min=max(10, 1e-4 * abs(prev_bic))


    return (curr_bic - prev_bic) <= -delta_bic_min

def ok_new_component(gmm, n, min_frac=0.03):
    w = gmm.weights_
    return (w.min() >= min_frac) and ((w * n).min() >= max(30, min_frac*n))



def determine_gmm_anchors(X_core,min_k=1,max_k=3,reg_covar=1e-6, random_state=0):
    """
    X_core: (n_core, 2) UMAP coords for core points
    Returns:
      labels: (n_core,) hard cluster labels in [0..k-1]
      centers: (k, 2) component means (anchors)
      model: fitted GaussianMixture
    """

    bic_list=[]

    n = len(X_core)
    if n == 0:
        raise ValueError("No core points provided.")
    max_k = max(min_k, min(max_k, n))  # cap by data size

    # Try k=1..max_k, pick best by BIC
    best_gmm = None
    k_pos=-1
    best_bic = np.inf
    prev_bic = None

    for k in range(min_k, max_k):
        gm = GaussianMixture(
            n_components=k,
            covariance_type="diag",      # 'diag' is even faster; 'full' is fine in 2D
            init_params="kmeans",
            n_init=3,
            max_iter=200,
            reg_covar=reg_covar,
            random_state=random_state,
        ).fit(X_core)
        bic = gm.bic(X_core)
        bic_list.append(bic)


        if k == 1 or bic < best_bic:
            # accept as current best
            best_bic, best_gmm = bic, gm

        if k > 1:
            if not should_continue(prev_bic, bic, delta_bic_min=10.0):
                break
            if not ok_new_component(gm, len(X_core), min_frac=0.03):
                break

        prev_bic = bic



    labels  = best_gmm.predict(X_core)
    centers = best_gmm.means_          # use as skeleton anchor points
    probs=best_gmm.predict_proba(X_core)

    return centers,labels,probs


#The resolution can be made parameter free in the visualization step.
def find_anchors(core_obj,anchor_finding_mode='binary',auto_select_resolution=False):


    X=core_obj.X
    n=X.shape[0]


    adj_list=core_obj.adj_list
    edges = [(u, v) for u in range(n) for v in adj_list[u] if u != v]

    G = ig.Graph(edges=edges, directed=True)
    G.vs["id"] = list(range(n))


    layers_=core_obj.layers_
    #densify=core_obj.densify
    densify=False

    core_nodes=layers_[0]

    G_core = G.induced_subgraph(core_nodes)

    # if densify:
    #     print("densifying...")
    #     G_core = densify_v0(adj_list, G, G_core, core_nodes)

    #Obtain the suitable resolution.
    if auto_select_resolution:
        res=choose_stopping_res(G,core_nodes,thr=0.95)

    else:
        res=core_obj.resolution

    partition = leidenalg.find_partition(
        G_core,
        leidenalg.RBConfigurationVertexPartition,  # modularity-based partition
        resolution_parameter=res,  # <-- control resolution here
        seed=np.random.randint(int(1e9))
    )



    filtered_partition=[]
    for cluster in partition:
        if len(cluster)>5:
            filtered_partition.append(cluster)

    partition=filtered_partition

    print("Total number of clusters=",len(partition))
    print([len(cluster) for cluster in partition])

    t1=time.time()



    anchors_dict = []
    for cluster in tqdm(partition, desc="Fitting GMM anchors"):
        anchors_dict.append(
            determine_gmm_anchors(X[np.array(list(cluster)).astype(int)])
        )

    #
    # anchors_dict = Parallel(n_jobs=6,backend='loky', batch_size='auto',prefer='processes')(
    #     delayed(determine_gmm_anchors)(X[np.array(list(partition[x])).astype(int)]) for x in tqdm(range(len(partition))))

    t2=time.time()


    print(f"GMM time={t2-t1:.3f} seconds (corrected)")

    anchor_list = []
    anchors = []
    t=0

    #Get the overall prob_matrix?


    if anchor_finding_mode=='default':

        for i,cluster in enumerate(partition):

            center_num=len(anchors_dict[i][0])

            ng_list=[[] for _ in range(center_num)]
            temp_labels=anchors_dict[i][1]
            n_1=len(temp_labels)
            for node in range(n_1):
                ng_list[temp_labels[node]].append(G_core.vs["id"][cluster[node]])




            for ng in ng_list:
                anchor_list.append(ng)
                anchors.append(np.mean(X[np.array(ng).astype(int)],axis=0))

        #Send the points connected to each center in the GMM-on-Leiden setup.

        final_prob_vec=None

    elif anchor_finding_mode=='smoothed':

        final_prob_vec = [{} for _ in range(n)]

        cx = 0
        for i,cluster in enumerate(partition):



            center_num=len(anchors_dict[i][0])
            ng_list=[[] for _ in range(center_num)]
            temp_labels=anchors_dict[i][1]
            n_1=len(temp_labels)
            for node in range(n_1):


                for i1 in range(center_num):
                    ng_list[i1].append(G_core.vs["id"][cluster[node]]) #Add the neighbors to each center.
                    final_prob_vec[G_core.vs["id"][cluster[node]]][cx+i1]=anchors_dict[i][2][node,i1] #For each node, add the tuples of center_number,probability. This will be used later. in create_anchor_edge_list

            cx+=center_num


            for ng in ng_list:
                anchor_list.append(ng)
                anchors.append(np.mean(X[np.array(ng).astype(int)],axis=0))


        tc=0



        print(len(anchor_list))
        for  u in range(len(anchor_list)):
            for v in anchor_list[u]:
                tc+=len(final_prob_vec[v].keys())

        print("Total GMM nodes=",tc)


    anchor_distances = []



    for i in range(len(anchor_list)):
        vec=[]
        for j in anchor_list[i]:
            vec.append(np.linalg.norm(anchors[i]-X[j]))

        anchor_distances.append(vec)

    anchored_list_sorted=[]
    anchored_dist_sorted=[]
    for idxs, dists in zip(anchor_list, anchor_distances):
        # zip indexes with distances, sort by distance, unzip back
        pairs = sorted(zip(idxs, dists), key=lambda x: x[1])
        idx_sorted, dist_sorted = zip(*pairs)

        anchored_list_sorted.append(list(idx_sorted))
        anchored_dist_sorted.append(list(dist_sorted))

    return anchored_list_sorted,anchored_dist_sorted,final_prob_vec





