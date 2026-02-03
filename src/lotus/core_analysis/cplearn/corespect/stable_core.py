from .stable_core_utils import find_intersection_layers
from ..utils.gen_utils import cutoffs

from ..utils.densify import densify_knn, densify_rw

import numpy as np
def stable_core(X,knn_list, flowrank_score, auto_select_core_frac,core_frac,ng_num,resolution,densification,**kwargs):

    n=len(knn_list)



    sorted_nodes= sorted(flowrank_score, key=lambda k: flowrank_score[k], reverse=True)

    if auto_select_core_frac:
        values = sorted(flowrank_score.values(), reverse=True)
        cuts = cutoffs(values)
        core_nodes = np.array(sorted_nodes[:int(min(cuts["cumu_mass"], 0.25 * n))]).astype(int)

    else:
        core_nodes=np.array(sorted_nodes[:int(core_frac*n)])

    if densification is not False:

        if densification=="rw":
            adj_list=densify_rw(knn_list, core_nodes)
        elif densification=="k-nn":
            adj_list=densify_knn(X,knn_list, core_nodes)
        else:
            raise ValueError(f"Invalid densification mode: {densification}")


        #Re-order the edges with

    else:
        adj_list=knn_list.copy()

    layer_candidate=find_intersection_layers(adj_list,core_nodes,ng_num=ng_num,rem_nodes=None,resolution=resolution)
    core_nodes_candidate = find_intersection_layers(adj_list, core_nodes, ng_num=ng_num, rem_nodes=np.array(core_nodes).astype(int), resolution=resolution)
    layer_candidate = list(core_nodes_candidate) + list(layer_candidate)



    return [layer_candidate],adj_list

