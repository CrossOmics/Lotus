#Main class for visualization.
import umap

from .viz_tools import anchored_map

import numpy as np

class Coremap:
    def __init__(self,
        corespect_obj,
        global_umap=None,
        fast_view=True,
        anchor_finding_mode='default',
        anchor_reach=None
    ):

        self.fast_view=fast_view
        self.anchor_finding_mode=anchor_finding_mode
        self.anchor_reach=anchor_reach


        self.core_obj=corespect_obj
        self.global_umap=global_umap
        self.X=corespect_obj.X
        self.X_umap=global_umap
        self.layers_=corespect_obj.layers_

        if self.global_umap is None:
            reducer = umap.UMAP(init='spectral')
            self.X_umap = reducer.fit_transform(self.X)


        if self.fast_view is True:
            label_dict = {}
            curr_layer = []
            for rounds, layer in enumerate(self.layers_):

                curr_layer.extend(layer)
                label_dict[rounds] = self.X_umap[np.array(curr_layer).astype(int)]

            self.label_dict = label_dict

        else:
            label_dict=anchored_map(self)


        self.label_dict=label_dict








