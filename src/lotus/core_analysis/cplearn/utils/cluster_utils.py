import numpy as np
import leidenalg

#TODO: Add choose_stopping_res function and add to config options
def make_subgraph(adj_list, core_nodes):
    import igraph as ig

    n = len(adj_list)
    edges = [(u, v) for u in range(n) for v in adj_list[u] if u != v]

    G = ig.Graph(edges=edges, directed=True)
    G.vs["id"] = list(range(n))  # Store original node ids

    G_core = G.induced_subgraph(core_nodes)

    return G,G_core


def cluster_subset(adj_list, core_nodes, resolution, **kwargs):

    """
    Cluster the core subgraph using the Leiden algorithm (directed version).

    Parameters
    ----------
    adj_list : list[list[int]]
        Adjacency list of the full directed graph.
    core_nodes : list[int]
        Indices of nodes forming the core subset.
    resolution : float
        Resolution parameter for Leiden clustering.
    **kwargs :
        Additional keyword arguments passed to leidenalg.find_partition.

    Returns
    -------
    labels : np.ndarray of shape (n,)
        Cluster labels for all nodes, with -1 for non-core nodes.
    """


    n=len(adj_list)
    G,G_core= make_subgraph(adj_list,core_nodes)

    partition = leidenalg.find_partition(
        G_core,
        leidenalg.RBConfigurationVertexPartition,  # modularity-based partition
        resolution_parameter=resolution,  # <-- control resolution here
        seed=np.random.randint(int(1e9)), **kwargs
    )

    id_map = np.array(G_core.vs["id"], dtype=int)

    labels = -1 * np.ones(n, dtype=int)
    labels[id_map] = partition.membership

    return labels