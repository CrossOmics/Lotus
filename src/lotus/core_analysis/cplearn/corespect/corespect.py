from dataclasses import asdict
from dataclasses import replace

from .config import (
    FlowRankConfig, StableCoreConfig,
    FineGrainedConfig, ClusterConfig, PropagationConfig
)


from .ranking import FlowRank
from .stable_core import stable_core
from .fine_grained_core import fine_grained_core
from .cluster_core import cluster_core
from .propagate import propagate_from_core


from ..utils.gen_utils import get_kNN

class CorespectModel:
    """
    CoreSpect model orchestrator.

    Stages:
        1. flowrank()              â node ranking and density estimation
        2. stable_core()           â stable core extraction
        3. fine_grained_core()     â optional refinement
        4. propagation_from_core() â diffusion or label spreading
    """

    def __init__(self, X,
                 flowrank_cfg=None,
                 stable_cfg=None,
                 fine_core_cfg=None,
                 cluster_cfg=None,
                 prop_cfg=None):


        self.X = X
        self.count_mat=None


        self.flowrank_cfg = flowrank_cfg or FlowRankConfig()
        self.stable_cfg = stable_cfg or StableCoreConfig()
        self.fine_core_cfg = fine_core_cfg or FineGrainedConfig()
        self.cluster_cfg = cluster_cfg or ClusterConfig()
        self.prop_cfg = prop_cfg or PropagationConfig()


        self.knn_list = None
        self.knn_dist = None


        self.adj_list = None
        self.flowrank_score_ = None
        self.core_layers_ = None
        self.fg_layers_= None
        self.prop_layers_ = None


        self.layers_ = None
        self.labels_ = None
        self.prob_matrix_ = None


        self.resolution=None

    def update_config(self, stage: str, **kwargs):
        """
        Update configuration parameters for a given stage.
        Example:
            model.update_config("flowrank", q=50, r=20)
        """
        valid = {
            "flowrank": "flowrank_cfg",
            "stable": "stable_cfg",
            "fine": "fine_core_cfg",
            "cluster": "cluster_cfg",
            "propagation": "prop_cfg",
        }
        if stage not in valid:
            raise ValueError(f"Unknown stage: {stage}. Valid: {list(valid)}")

        cfg_attr = valid[stage]
        old_cfg = getattr(self, cfg_attr)
        new_cfg = replace(old_cfg, **kwargs)
        setattr(self, cfg_attr, new_cfg)
        print(f" Updated {stage} config: {kwargs}")

        return self




    def build_graph(self):
        """Construct adjacency list based on the largest neighborhood size."""
        q = self.flowrank_cfg.q
        r = self.flowrank_cfg.r
        ng = self.stable_cfg.ng_num
        k = max(q, r, ng)

        knn_list,knn_dist = get_kNN(self.X,k)

        self.knn_list = knn_list
        self.knn_dist = knn_dist

        return self

    def flowrank(self):
        self.flowrank_score_ = FlowRank(self.knn_list, **asdict(self.flowrank_cfg))
        return self

    def stable_core(self):
        self.core_layers_,self.adj_list = stable_core(self.X,self.knn_list, self.flowrank_score_, **asdict(self.stable_cfg))
        self.layers_ = self.core_layers_

        self.resolution = self.stable_cfg.resolution

        return self


    def fine_grained_core(self):
        self.fg_layers_ =fine_grained_core(self.X, self.adj_list, self.core_layers_, **asdict(self.fine_core_cfg), **asdict(self.flowrank_cfg))
        self.layers_=self.fg_layers_

        return self

    def cluster_core(self):


        if self.fg_layers_ is None:
            core_nodes = self.core_layers_[0]
        else:
            core_nodes = self.fg_layers_[0]

        self.labels_ =cluster_core(self.X,self.adj_list,core_nodes,**asdict(self.cluster_cfg))


    def propagation_from_core(self):

        if self.fg_layers_ is None:
            layers = self.core_layers_

        else:
            print([len(layer) for layer in self.fg_layers_],flush=True)
            layers = self.fg_layers_

        self.labels_,self.prop_layers_,self.prob_matrix_ = propagate_from_core(self.adj_list,layers,self.labels_
        , **asdict(self.prop_cfg)
        )

        self.layers_=self.prop_layers_


        return self

    def run(self, fine_grained=False, propagate=True):
        self.build_graph()
        self.flowrank()
        self.stable_core()

        if fine_grained:

            print("Finding Fine Grained Core",flush=True)
            self.fine_grained_core()

        self.cluster_core()

        if propagate:
             self.propagation_from_core()

        return self

