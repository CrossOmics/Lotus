import numpy as np
import plotly.graph_objects as go
import plotly.express as px

def visualize_coremap(coremap_obj, labels=None, use_webgl=True):
    """
    Plot layered 2D embeddings with a discrete slider.
    - One trace per cluster per layer.
    - Legend updates dynamically for the active layer only.
    """

    layers = []
    for coord_local in coremap_obj.label_dict:
        layers.append(coremap_obj.label_dict[coord_local])

    layer_indices = []
    curr_layer_local = []
    for layer_local in coremap_obj.layers_:
        curr_layer_local.extend(layer_local)
        layer_indices.append(curr_layer_local.copy())


    if labels is None:
        n_t=coremap_obj.X.shape[0]
        labels=np.zeros(n_t)

    else:
        labels = np.asarray(labels)

    layer_indices = [np.asarray(idxs) for idxs in layer_indices]

    # ---- color map (stable across layers) ----
    unique_labels = np.unique(labels)
    palette = (
        px.colors.qualitative.Plotly
        + px.colors.qualitative.D3
        + px.colors.qualitative.Set1
        + px.colors.qualitative.Set2
        + px.colors.qualitative.Set3
        + px.colors.qualitative.T10
        + px.colors.qualitative.Dark24
        + px.colors.qualitative.Light24
    )
    color_map = {lab: palette[i % len(palette)] for i, lab in enumerate(sorted(unique_labels))}

    # ---- axis ranges (fixed across all layers) ----
    all_x = np.concatenate([emb[:, 0] for emb in layers])
    all_y = np.concatenate([emb[:, 1] for emb in layers])
    pad_x = 0.05 * (all_x.max() - all_x.min())
    pad_y = 0.05 * (all_y.max() - all_y.min())
    x_range = [all_x.min() - pad_x, all_x.max() + pad_x]
    y_range = [all_y.min() - pad_y, all_y.max() + pad_y]

    traces = []
    layer_trace_indices = []

    # ---- build traces: one per cluster per layer ----
    for layer_id, (emb, idxs) in enumerate(zip(layers, layer_indices)):
        start = len(traces)

        for lab in unique_labels:
            mask = labels[idxs] == lab
            if not np.any(mask):
                continue

            if use_webgl:
                Trace = go.Scattergl
                marker_style = dict(
                    color=color_map[lab],
                    size=6,
                    opacity=0.7
                )
            else:
                Trace = go.Scatter
                marker_style = dict(
                    color=color_map[lab],
                    size=9,
                    opacity=0.85,
                    line=dict(width=1, color="white"),  # halo effect
                    symbol="circle"
                )

            traces.append(
                Trace(
                    x=emb[mask, 0],
                    y=emb[mask, 1],
                    mode="markers",
                    marker=marker_style,
                    name=str(lab),
                    showlegend=(layer_id == 0),   # only show legend for layer 0 initially
                    visible=(layer_id == 0)       # only layer 0 visible initially
                )
            )

        end = len(traces)
        layer_trace_indices.append((start, end))

    # ---- slider steps ----
    steps = []
    n_traces = len(traces)
    for k, (start, end) in enumerate(layer_trace_indices):
        vis = [False] * n_traces
        showleg = [False] * n_traces

        for i in range(start, end):
            vis[i] = True
            showleg[i] = True

        steps.append(dict(
            method="update",
            args=[{"visible": vis, "showlegend": showleg},
                  {"title": f"Layer {k}"}],
            label=str(k)
        ))

    sliders = [dict(active=0, steps=steps)]

    # ---- final figure ----
    fig = go.Figure(
        data=traces,
        layout=go.Layout(
            template="plotly_white",
            width=1000, height=800,
            xaxis=dict(range=x_range, title="Dim 1"),
            yaxis=dict(range=y_range, title="Dim 2"),
            sliders=sliders,
            title="Embedding Layers with Dynamic Legend"
        )
    )
    return fig
