# produce_umap_plots.py
#
# This file produces UMAP plots (figures 3, 4, and 5 in the manuscript)
#
#
# This file requires the umap module:
#     https://umap-learn.readthedocs.io/en/latest/index.html
#     https://github.com/lmcinnes/umap
# It also requires plotly:
#     https://plotly.com/python/getting-started/
#
# This file takes less than a minute to run.
#
# Input files:
#     ./results38/hg38_chr6_200datacorrelation.h5
# (Correlation distance for chromosome of interest)
#
# To change the chromosome, change the file path for the correlation matrix
# and change the chr_id_1 variable.
#
# The script uses plotly to create an interactive html file which
# can be opened in an ordinary web browser.
#
# Output files:
#     cells38_nolegend.html
#     epi38_nolegend.html
#     epi38_histonevsnot.html
#     epi38_activity_sh.html
#     epi38_factor_sh.html
#     epi38_activity_sh_legend.png
#     epi38_factor_sh_legend.png
#
# For a description of the content of these plots, see figures 3, 4 and 5 in the manuscript.
#
import pandas as pd
import numpy as np

# Note: The following allows you to interact with the
# plotly plots in a jupyter notebook.
from IPython.core.interactiveshell import InteractiveShell

InteractiveShell.ast_node_interactivity = "all"

from ast import literal_eval
import plotly.offline as py

py.init_notebook_mode(connected=True)
import plotly.graph_objs as go
from umap import UMAP

import random

marks = pd.read_excel("./data/target_activity_factor.xlsx")
marks.columns = ["Target", "Activity", "Factor"]
marks.head(5)
marks["Activity"].value_counts()
marks["Factor"].value_counts()

df38 = pd.read_csv("./data/genome_df38.csv", delimiter=",")
df38 = pd.DataFrame(df38)
df38 = df38.loc[:, ["Accession", "Target", "Biosample term name", "Genome"]]
df38["Target"].value_counts()
df38["Biosample term name"].value_counts()
merged_df = pd.merge(df38, marks, on="Target")
merged_df.head(3)
df38 = merged_df
df38["class"] = df38["Target"].apply(lambda x: "Histone" if x.startswith("H") and x != "HDAC2" else "Non-Histone")
df38["activity_class"] = df38["class"] + str(": ") + df38["Activity"]
df38["factor_class"] = df38["class"] + str(": ") + df38["Factor"]

# Figure 3
chr_id_1 = 6
cor_dist = pd.read_csv(
    "./results38/hg38_chr" + str(chr_id_1) + "_200data" + "correlation.h5", index_col=0
)  # pre-computed cor for chr 6 hg38 (look at the script that producing it)
cor_dist = cor_dist.to_numpy()
embedding = UMAP(n_neighbors=10, n_components=3, metric="precomputed", random_state=42).fit_transform(cor_dist)
dfu = pd.DataFrame(embedding, columns=("x", "y", "z"))
dfu.index = df38["Biosample term name"]
dfu["class"] = dfu.index
labels = np.unique(dfu["class"])

colors = [
    "rgb(0, 0, 128)",
    "rgb(128, 128, 128)",
    "rgb(0, 0, 0)",
    "rgb(230, 25, 75)",
    "rgb(60, 180, 75)",
    "rgb(255, 225, 25)",
    "rgb(0, 130, 200)",
    "rgb(245, 130, 48)",
    "rgb(145, 30, 180)",
    "rgb(70, 240, 240)",
    "rgb(240, 50, 230)",
    "rgb(210, 245, 60)",
    "rgb(250, 190, 212)",
    "rgb(0, 128, 128)",
    "rgb(220, 190, 255)",
    "rgb(170, 110, 40)",
    "rgb(255, 250, 200)",
    "rgb(128, 0, 0)",
    "rgb(170, 255, 195)",
    "rgb(128, 128, 0)",
    "rgb(255, 215, 180)",
    "rgb(255, 255, 255)",
]
data_graph = []
for no, name in enumerate(np.unique(dfu["class"])):
    graph = go.Scatter3d(
        x=dfu[dfu["class"] == name]["x"],
        y=dfu[dfu["class"] == name]["y"],
        z=dfu[dfu["class"] == name]["z"],
        name=labels[no],
        mode="markers",
        marker=dict(
            size=8,
            color="#%02x%02x%02x" % literal_eval(colors[no][3:]),  #'#%02x%02x%02x' % literal_eval(colors[no][3:]),
            line=dict(color=literal_eval(colors[no][3:]), width=0.5),  #            color = '#%02x%02x%02x' % literal_eval(colors[no][3:]),
            opacity=1,
        ),
    )
    data_graph.append(graph)

layout = go.Layout(scene=dict(camera=dict(eye=dict(x=0.5, y=0.5, z=0.5))), margin=dict(l=0, r=0, b=0, t=0), showlegend=False)  # hide legend
fig = go.Figure(data=data_graph, layout=layout)
py.iplot(fig, filename="3d-scatter")
fig.write_html("cells38_nolegend.html")


# Figure 4
dfu.index = df38["Target"]
dfu["class"] = dfu.index
labels = np.unique(dfu["class"])


colors = [
    "rgb(230, 25, 75)",
    "rgb(60, 180, 75)",
    "rgb(255, 225, 25)",
    "rgb(0, 130, 200)",
    "rgb(245, 130, 48)",
    "rgb(145, 30, 180)",
    "rgb(70, 240, 240)",
    "rgb(240, 50, 230)",
    "rgb(210, 245, 60)",
    "rgb(250, 190, 212)",
    "rgb(0, 128, 128)",
    "rgb(220, 190, 255)",
    "rgb(170, 110, 40)",
    "rgb(255, 250, 200)",
    "rgb(128, 0, 0)",
    "rgb(170, 255, 195)",
    "rgb(128, 128, 0)",
    "rgb(255, 215, 180)",
    "rgb(0, 0, 128)",
    "rgb(128, 128, 128)",
    "rgb(255, 255, 255)",
    "rgb(0, 0, 0)",
    "rgb(0,0,0)",
    "rgb(255,0,0)",
    "rgb(0,255,0)",
    "rgb(0,0,255)",
    "rgb(255,255,0)",
    "rgb(0,255,255)",
    "rgb(255,0,255)",
    "rgb(192,192,192)",
    "rgb(128,128,128)",
    "rgb(128,0,0)",
    "rgb(128,128,0)",
    "rgb(0,128,0)",
    "rgb(128,0,128)",
    "rgb(0,128,128)",
    "rgb(0,0,128)",
    "rgb(250,128,114)",
    "rgb(255,160,122)",
    "rgb(255,69,0)",
    "rgb(255,140,0)",
    "rgb(255,165,0)",
    "rgb(255,215,0)",
    "rgb(184,134,11)",
    "rgb(218,165,32)",
    "rgb(238,232,170)",
    "rgb(189,183,107)",
    "rgb(240,230,140)",
    "rgb(128,128,0)",
    "rgb(255,255,0)",
    "rgb(154,205,50)",
    "rgb(85,107,47)",
    "rgb(107,142,35)",
    "rgb(124,252,0)",
    "rgb(34,139,34)",
    "rgb(0,100,0)",
    "rgb(32,178,170)",
    "rgb(47,79,79)",
    "rgb(0,128,128)",
    "rgb(176,224,230)",
    "rgb(95,158,160)",
    "rgb(70,130,180)",
    "rgb(100,149,237)",
    "rgb(0,191,255)",
    "rgb(30,144,255)",
    "rgb(0,0,139)",
    "rgb(0,0,205)",
    "rgb(65,105,225)",
    "rgb(138,43,226)",
    "rgb(75,0,130)",
    "rgb(72,61,139)",
    "rgb(106,90,205)",
    "rgb(123,104,238)",
    "rgb(147,112,219)",
    "rgb(139,0,139)",
    "rgb(148,0,211)",
    "rgb(255,0,255)",
    "rgb(218,112,214)",
    "rgb(199,21,133)",
    "rgb(219,112,147)",
    "rgb(255,20,147)",
    "rgb(255,105,180)",
    "rgb(255,192,203)",
    "rgb(250,235,215)",
    "rgb(245,245,220)",
    "rgb(255,228,196)",
    "rgb(160,82,45)",
    "rgb(210,105,30)",
    "rgb(205,133,63)",
    "rgb(244,164,96)",
    "rgb(222,184,135)",
    "rgb(210,180,140)",
    "rgb(188,143,143)",
    "rgb(245,255,250)",
    "rgb(112,128,144)",
    "rgb(176,196,222)",
    "rgb(230,230,250)",
]
colors = np.unique(colors)
colors = random.choices(colors, k=len(np.unique(dfu.index)))

data_graph = []
for no, name in enumerate(np.unique(dfu["class"])):
    graph = go.Scatter3d(
        x=dfu[dfu["class"] == name]["x"],
        y=dfu[dfu["class"] == name]["y"],
        z=dfu[dfu["class"] == name]["z"],
        name=labels[no],
        mode="markers",
        marker=dict(
            size=8,
            color="#%02x%02x%02x" % literal_eval(colors[no][3:]),  #'#%02x%02x%02x' % literal_eval(colors[no][3:]),
            line=dict(color=literal_eval(colors[no][3:]), width=0.5),  #            color = '#%02x%02x%02x' % literal_eval(colors[no][3:]),
            opacity=1,
        ),
    )
    data_graph.append(graph)

layout = go.Layout(scene=dict(camera=dict(eye=dict(x=0.5, y=0.5, z=0.5))), margin=dict(l=0, r=0, b=0, t=0), showlegend=False)  # hide legend
fig = go.Figure(data=data_graph, layout=layout)
py.iplot(fig, filename="3d-scatter")
fig.write_html("epi38_nolegend.html")

# Histone vs non-Histone UMAP plot
import pandas as pd
import plotly.graph_objs as go

# camera = dict(
#     eye=dict(x=1.5, y=1.5, z=1.5)
# )
dfu = pd.DataFrame(embedding, columns=("x", "y", "z"))
dfu.index = df38["class"]
# dfu["class"] = df38['Type']
# Assign colors: Red for Histone and Blue for Non-Histone
color_map = {"Histone": "blue", "Non-Histone": "green"}
dfu["color"] = dfu.index.map(color_map)

# Plotting
data_graph = []
for class_name, color in color_map.items():
    mask = dfu.index == class_name
    graph = go.Scatter3d(
        x=dfu[mask]["x"],
        y=dfu[mask]["y"],
        z=dfu[mask]["z"],
        name=class_name,
        mode="markers",
        marker=dict(size=8, color=color, line=dict(color=color, width=0.5), opacity=1),
    )
    data_graph.append(graph)

layout = go.Layout(scene=dict(camera=dict(eye=dict(x=0.5, y=0.5, z=0.5))), margin=dict(l=0, r=0, b=0, t=0))
fig = go.Figure(data=data_graph, layout=layout)
fig.show()
fig.write_html("epi38_histonevsnot.html")

# Figure 5 Upper
dfu.index = df38["activity_class"]
color_map = {
    "Non-Histone: Permissive": "red",
    "Histone: Permissive": "blue",
    "Non-Histone: Repressive": "green",
    "Non-Histone: Both": "orange",
    "Histone: Repressive": "purple",
    "Non-Histone: Both*": "pink",
    "Histone: H2A histone gene": "fuchsia",
    "Histone: Histone 3 gene": "grey",
}
dfu["color"] = dfu.index.map(color_map)

marker_shape_map = {"Histone": "cross", "Non-Histone": "circle"}


def get_marker_shape(type_value):
    return marker_shape_map[type_value]


dfu.index = df38["class"]
dfu["marker_shape"] = dfu.index.map(get_marker_shape)


data_graph = []
for type_value in dfu.index.unique():
    mask = dfu.index == type_value
    graph = go.Scatter3d(
        x=dfu[mask]["x"],
        y=dfu[mask]["y"],
        z=dfu[mask]["z"],
        name=type_value,
        mode="markers",
        marker=dict(
            size=8, color=dfu[mask]["color"], symbol=dfu[mask]["marker_shape"], line=dict(color=dfu[mask]["color"], width=0.5), opacity=1
        ),
    )
    data_graph.append(graph)

fig = go.Figure(data=data_graph, layout=layout)
fig.show()
fig.write_html("epi38_activity_sh.html")

# Figure 5 Upper legend
import plotly.graph_objs as go
import plotly.io as pio

categories = {
    "Non-Histone: Permissive": {"color": "red", "symbol": "circle"},
    "Non-Histone: Repressive": {"color": "green", "symbol": "circle"},
    "Non-Histone: Both": {"color": "orange", "symbol": "circle"},
    "Non-Histone: Both*": {"color": "pink", "symbol": "circle"},
    "Histone: Permissive": {"color": "blue", "symbol": "cross"},
    "Histone: Repressive": {"color": "purple", "symbol": "cross"},
    "Histone: H2A histone gene": {"color": "fuchsia", "symbol": "cross"},
    "Histone: Histone 3 gene": {"color": "grey", "symbol": "cross"},
}


traces = [
    go.Scatter(
        x=[0], y=[0], mode="markers", marker=dict(size=10, color=props["color"], symbol=props["symbol"]), name=category  # Dummy data
    )
    for category, props in categories.items()
]

fig = go.Figure(data=traces)

fig.update_layout(
    margin=dict(l=0, r=0, t=0, b=0),  # Minimal margins
    xaxis=dict(showgrid=False, showticklabels=False, zeroline=False, visible=False),
    yaxis=dict(showgrid=False, showticklabels=False, zeroline=False, visible=False),
    showlegend=True,
    legend=dict(orientation="h", x=0.5, y=0.5, xanchor="center", yanchor="middle", title="Activity:"),
    plot_bgcolor="rgba(0,0,0,0)",
)

file_path = "epi38_activity_sh_legend.png"
fig.write_image(file_path, format="png")

# Figure 5 Lower Plot
import plotly.graph_objs as go
import plotly.io as pio
import pandas as pd
import plotly.graph_objs as go

df38["factor_class"] = df38["factor_class"].replace(
    "Non-Histone: post-transcriptional regulator", "Non-Histone: Post-transcriptional regulator"
)

dfu = pd.DataFrame(embedding, columns=("x", "y", "z"))
dfu.index = df38["factor_class"]
color_map = {
    "Non-Histone: Transcription factor": "crimson",
    "Non-Histone: Chromatin modifier": "darkblue",
    "Non-Histone: Histone modifier": "lightgreen",
    "Non-Histone: Post-transcriptional regulator": "magenta",
    "Histone: Histone modification": "gold",
    "Histone: Histone protein": "plum",
}
dfu["color"] = dfu.index.map(color_map)

marker_shape_map = {"Histone": "cross", "Non-Histone": "circle"}


def get_marker_shape(type_value):
    return marker_shape_map[type_value]


dfu.index = df38["class"]
dfu["marker_shape"] = dfu.index.map(get_marker_shape)

data_graph = []
for type_value in dfu.index.unique():
    mask = dfu.index == type_value
    graph = go.Scatter3d(
        x=dfu[mask]["x"],
        y=dfu[mask]["y"],
        z=dfu[mask]["z"],
        name=type_value,
        mode="markers",
        marker=dict(
            size=8, color=dfu[mask]["color"], symbol=dfu[mask]["marker_shape"], line=dict(color=dfu[mask]["color"], width=0.5), opacity=1
        ),
    )
    data_graph.append(graph)

fig = go.Figure(data=data_graph, layout=layout)
fig.show()
fig.write_html("epi38_factor_sh.html")

# Figure 5 Lower legend

categories = {
    "Non-Histone: Transcription factor": {"color": "crimson", "symbol": "circle"},
    "Non-Histone: Chromatin modifier": {"color": "darkblue", "symbol": "circle"},
    "Non-Histone: Histone modifier": {"color": "lightgreen", "symbol": "circle"},
    "Non-Histone: Post-transcriptional regulator": {"color": "magenta", "symbol": "circle"},
    "Histone: Histone modification": {"color": "gold", "symbol": "cross"},
    "Histone: Histone protein": {"color": "plum", "symbol": "cross"},
}


traces = [
    go.Scatter(x=[0], y=[0], mode="markers", marker=dict(size=10, color=props["color"], symbol=props["symbol"]), name=category)
    for category, props in categories.items()
]

fig = go.Figure(data=traces)

fig.update_layout(
    margin=dict(l=0, r=0, t=0, b=0),  # Minimal margins
    xaxis=dict(showgrid=False, showticklabels=False, zeroline=False, visible=False),
    yaxis=dict(showgrid=False, showticklabels=False, zeroline=False, visible=False),
    showlegend=True,
    legend=dict(orientation="h", x=0.5, y=0.5, xanchor="center", yanchor="middle", title="Factor:"),
    plot_bgcolor="rgba(0,0,0,0)",
)

file_path = "epi38_factor_sh_legend.png"
fig.write_image(file_path, format="png")

file_path
