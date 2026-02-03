from .stable_core_utils import choose_stopping_res, cluster_subset
from ..utils.densify import densify_knn, densify_rw





def cluster_core(X,adj_list,core_nodes,densification,resolution,auto_select_resolution):

    if densification is not False:

        if densification=="rw":
            adj_list=densify_rw(adj_list, core_nodes)
        elif densification=="k-nn":
            adj_list=densify_knn(X,adj_list, core_nodes)
        else:
            raise ValueError(f"Invalid densification mode: {densification}")



    if auto_select_resolution:
        resolution= choose_stopping_res(adj_list, core_nodes, starting_resolution=resolution)



    final_labels=cluster_subset(adj_list,core_nodes,resolution)


    return final_labels