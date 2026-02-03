
from ..utils.densify import densify_knn, densify_rw
from ..utils.gen_utils import cutoffs
from .stable_core_utils import find_intersection_layers

from .ranking import FlowRank

import numpy as np

def fine_grained_core(X,adj_list,layers_,stable_core_densification,auto_select_core_frac,r,core_frac,fine_grain_densification,starting_resolution,**kwargs):

    stable_core = layers_[0]





    n=len(adj_list)
    n_core=len(stable_core)

    if stable_core_densification is not False:

        print("Densifying the induced subgraph of stable core points")

        if stable_core_densification=="rw":
            adj_list_fg=densify_rw(adj_list, stable_core)
        elif stable_core_densification=="k-nn":
            adj_list_fg=densify_knn(X,adj_list, stable_core)
        else:
            raise ValueError(f"Invalid densification mode: {stable_core_densification}")

    else:
        adj_list_fg=adj_list.copy()



    # ---- ---- ---- ---- Everything below is limited to stable core points ---- ---- ---- ----#

    stable_map=-1*np.ones(n,dtype=int)
    for i,node in enumerate(stable_core):
        stable_map[node]=i



    adj_list_core=[[] for _ in range(n_core)]
    for u in stable_core:
        adj_list_core[stable_map[u]]=[stable_map[v] for v in adj_list_fg[u] if stable_map[v] != -1]



    X_core = X[stable_core]



    #--- Get FlowRank on the core points ----

    flowrank_score=FlowRank(adj_list_core,r)
    sorted_nodes= sorted(flowrank_score, key=lambda k: flowrank_score[k], reverse=True)

    if auto_select_core_frac:
        values = sorted(flowrank_score.values(), reverse=True)
        cuts = cutoffs(values)
        core_nodes = np.array(sorted_nodes[:int(min(cuts["cumu_mass"], 0.25 * n))]).astype(int)

    else:
        core_nodes=np.array(sorted_nodes[:int(core_frac*n)])


    #---- Densify the core points again ----#

    if fine_grain_densification is not False:

        print("Densifying the induced subgraph of fine grained core points")

        if fine_grain_densification=="rw":
            adj_list_fgc=densify_rw(adj_list_core, core_nodes)
        elif fine_grain_densification=="k-nn":
            adj_list_fgc=densify_knn(X_core,adj_list_core, core_nodes)
        else:
            raise ValueError(f"Invalid densification mode: {fine_grain_densification}")

    else:
        adj_list_fgc=adj_list_core.copy()


    #--- We can select initial resolution here with our filtering strategy ----#
    #For now let's just go with starting resolution

    #Fake ng_num. Currently, not being used in any case. Set to r



    layer_candidate=find_intersection_layers(adj_list_fgc,core_nodes,ng_num=r,rem_nodes=np.array(core_nodes).astype(int),resolution=starting_resolution)

    #This is the new core. Now we can update the layers.

    print("Number of original and filtered fine grained core nodes:",len(core_nodes),len(layer_candidate),flush=True)




    layer_core=[stable_core[u] for u in layer_candidate]
    layer_remaining = list(set(stable_core) - (set(layer_core)))

    new_layers_ = [layer_core, layer_remaining]

    print("Finished running fine_grained_core",flush=True)

    return new_layers_


