import networkit as nk
import networkx as nx
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from random import sample

G_nk = nk.readGraph("facebook_combined.txt", nk.Format.SNAP)
G = nk.nxadapter.nk2nx(G_nk)

# num_to_remove = int(len(G) / 3)
# nodes = sample(list(G.nodes), num_to_remove)
# G.remove_nodes_from(nodes)

# remove low-degree nodes
# low_degree = [n for n, d in G.degree() if d < 10]
# G.remove_nodes_from(low_degree)
# largest connected component

components = nx.connected_components(G)
largest_component = max(components, key=len)
H = G.subgraph(largest_component)

# compute centrality
centrality = dict(nx.betweenness_centrality(H, k=10, endpoints=True))
# compute community structure
lpc = nx.community.label_propagation_communities(H)
community_index = {n: i for i, com in enumerate(lpc) for n in com}


#### draw graph ####
fig, ax = plt.subplots(figsize=(20, 15),dpi=200)
# pos = nx.spring_layout(H, k=0.15, seed=4572321)
pos = nx.layout.kamada_kawai_layout(H)
node_color = [community_index[n] for n in H]
node_size = [v * 20000 for v in centrality.values()]


nx.draw_networkx_edges(
    H,
    pos=pos,
    edge_color='gainsboro',
    alpha=0.2,
)
nx.draw_networkx_nodes(
    H,
    pos=pos,
    node_color=node_color,
    node_size=node_size,
    alpha=0.4    
)

# Title/legend
# font = {"color": "k", "fontweight": "bold", "fontsize": 20}
# ax.set_title("Gene functional association network (C. elegans)", font)
# # Change font color for legend
# font["color"] = "r"

# ax.text(
#     0.80,
#     0.10,
#     "node color = community structure",
#     horizontalalignment="center",
#     transform=ax.transAxes,
#     fontdict=font,
# )
# ax.text(
#     0.80,
#     0.06,
#     "node size = betweeness centrality",
#     horizontalalignment="center",
#     transform=ax.transAxes,
#     fontdict=font,
# )

# Resize figure for label readibility
ax.margins(0.1, 0.05)
fig.tight_layout()
plt.axis("off")
# plt.savefig('fb_network_full_nolog.pdf')
# plt.savefig('fb_network_full_nolog.png')


