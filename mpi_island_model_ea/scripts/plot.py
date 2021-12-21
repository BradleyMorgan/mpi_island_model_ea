#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import os
import numpy as np
import importlib.util

ex_size='256'
ex_type='evo'

trunk = os.getcwd()
topos = trunk + '/results/logs/' + ex_size + '/' + ex_type + '/meta_001/stats/topos'

for entry in os.scandir(topos):

   if entry.is_file():

      print(entry.name)
      spec = importlib.util.spec_from_file_location(entry.name, entry.path)
      module = importlib.util.module_from_spec(spec)
      spec.loader.exec_module(module)
      G = nx.from_numpy_matrix(np.array(module.matrix), create_using=nx.MultiDiGraph())
      pos = nx.layout.spiral_layout(G)

      node_sizes = [3 + 10 * i for i in range(len(G))]
      M = G.number_of_edges()
      edge_colors = range(2, M + 2)
      edge_alphas = [(5 + i) / (M + 4) for i in range(M)]

      plt.figure(figsize=(16,12))

      nodes = nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color="blue")
      edges = nx.draw_networkx_edges(G, pos, node_size=node_sizes, arrowstyle="->", arrowsize=10, edge_color=edge_colors, edge_cmap=plt.cm.Blues, width=2)

      labels = {i : i + 1 for i in G.nodes()}
      nx.draw_networkx_labels(G, pos, labels, font_size=8,font_color='white')

      # set alpha value for each edge
      for i in range(M):
         edges[i].set_alpha(edge_alphas[i])

      pc = mpl.collections.PatchCollection(edges, cmap=plt.cm.Blues)
      pc.set_array(edge_colors)
      plt.colorbar(pc)

      ax = plt.gca()
      ax.set_axis_off()
      plt.show()

