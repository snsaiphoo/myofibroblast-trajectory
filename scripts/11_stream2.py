import scanpy as sc
import stream2 as st2
import pacmap
from pathlib import Path

FIG_DIR = Path("../results/stream2")
FIG_DIR.mkdir(exist_ok=True)

st2.settings.set_workdir(str(FIG_DIR))

adata = sc.read_h5ad("../data/mesenchymal_scored.h5ad")

pac = pacmap.PaCMAP(
    n_components=2,
    random_state=42
)

adata.obsm["X_pacmap"] = pac.fit_transform(
    adata.obsm["X_pca"]
)

adata.obsm["X_dr"] = adata.obsm["X_pacmap"]

print(adata.obsm.keys())

import matplotlib.pyplot as plt

sc.pl.embedding(
    adata,
    basis="X_pacmap",
    color="cell_type_refined",
    frameon=False,
    show=False
)

plt.savefig(
    FIG_DIR / "00_pacmap_celltypes.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# compute density weights 
st2.tl.get_weights(
    adata,
    obsm="X_dr",
    bandwidth=1,
    griddelta=100
)

# weighted seed graph
st2.tl.seed_graph(
    adata,
    obsm="X_dr",
    clustering="kmeans",
    n_clusters=30,
    use_weights=True
)

# learn the principal graph

st2.tl.learn_graph(
    adata,
    obsm="X_dr",
    n_nodes=50,
    use_seed=True,
    use_weights=True
)


# save the graph 
st2.pl.graph(
    adata,
    key="epg",
    color=["cell_type_refined"],
    show_node=True,
    fig_size=(10,8),
    save_fig=True,
    fig_path=str(FIG_DIR),
    fig_name="01_stream_graph_weighted.png"
)

# find the mature tenocyte node
st2.pl.graph(
    adata,
    key="epg",
    color=["cell_type_refined"],
    show_node=True,
    show_text=True,
    fig_size=(10,8),
    save_fig=True,
    fig_path=str(FIG_DIR),
    fig_name="02_graph_node_numbers.png"
)

# node 17 was chosen
ROOT_NODE = 17

st2.tl.infer_pseudotime(
    adata,
    source=ROOT_NODE
)

print([x for x in adata.obs.columns if "pseudo" in x.lower()])

# plot pseudotime
st2.pl.graph(
    adata,
    key="epg",
    color=["epg_pseudotime"],
    fig_size=(10,8),
    show_node=False,
    save_fig=True,
    fig_path=str(FIG_DIR),
    fig_name="03_stream_pseudotime.png"
)

# find the available plotting arguments
import inspect

print(inspect.signature(st2.pl.stream))

# stream plot colored by cell types
st2.pl.stream(
    adata,
    source=17,
    color=["cell_type_refined"],
    fig_size=(12,8),
    save_fig=True,
    fig_path=str(FIG_DIR),
    fig_format="png"
)

# then by pseudotime
st2.pl.stream(
    adata,
    source=17,
    color=["epg_pseudotime"],
    fig_size=(12,8),
    save_fig=True,
    fig_path=str(FIG_DIR),
    fig_format="png"
)

# Plot UCell Overlays 
ucell_scores = [
    "YAP_TEAD_UCell",
    "Integrin_FAK_UCell",
    "RhoA_ROCK_UCell",
    "Myofibroblast_UCell",
    "GOBP_FIBROBLAST_ACTIVATION_UCell",
    "GOBP_EXTRACELLULAR_STRUCTURE_ORGANIZATION_UCell",
    "HALLMARK_TGF_BETA_SIGNALING_UCell",
    "HALLMARK_INFLAMMATORY_RESPONSE_UCell"
]

for score in ucell_scores:
    print(f"Plotting {score}")
    st2.pl.graph(
        adata,
        key="epg",
        color=[score],
        fig_size=(10,8),
        show_node=False,
        save_fig=True,
        fig_path=str(FIG_DIR),
        fig_name=f"{score}.png"
    )

# plot the cell types one more time 
st2.pl.graph(
    adata,
    key="epg",
    color=["cell_type_refined"],
    fig_size=(10,8),
    show_node=False,
    save_fig=True,
    fig_path=str(FIG_DIR),
    fig_name="CellTypes.png"
)

st2.pl.graph(
    adata,
    key="epg",
    color=["epg_pseudotime"],
    fig_size=(10,8),
    show_node=False,
    save_fig=True,
    fig_path=str(FIG_DIR),
    fig_name="Pseudotime.png"
)
