#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 14:10:33 2023

@author: yutasuzuki
"""

import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sns


for i in np.arange(10):
    plt.figure()
    
    normal_data = np.random.randn(10, 10)
    # sns.heatmap(normal_data,cmap="coolwarm")
    # ax=sns.heatmap(normal_data,cmap="RdBu", annot=False)
    ax = sns.heatmap(normal_data,cmap="binary", annot=False)
    ax.set_aspect('equal')
    
    plt.savefig("img" + str(i) + ".pdf")

plt.figure()
for i in np.arange(6):
    
    normal_data = np.random.rand(30, 1)+i
    ax = plt.plot(np.arange(len(normal_data)),normal_data,'k')

plt.axis("off")
   
plt.savefig("line1.pdf")

#%%
import networkx as nx
import string
import itertools

G = nx.Graph()

layer = [5,6,7,7,6,5]

nodeName = list(string.ascii_lowercase) + [s+'0' for s in list(string.ascii_letters)]

nodeCount = 0
pos={}
for iLayer,l in enumerate(layer):
    
    for iNeuron in np.arange(l):
        G.add_node(nodeName[nodeCount])
        
        tmp_pos = np.linspace(0,l,l) + (max(layer) - l)/2
        
        pos[nodeName[nodeCount]] = [iLayer,tmp_pos[iNeuron]]

        nodeCount+=1


for iLayer,l in enumerate(layer[:-1]):
    l1 = nodeName[sum(layer[:iLayer]): sum(layer[:iLayer]) + layer[iLayer]]
    l2 = nodeName[sum(layer[:iLayer+1]): sum(layer[:iLayer+1]) + layer[iLayer+1]]
    for p in list(itertools.product(np.arange(layer[iLayer]), np.arange(layer[iLayer+1]))):
        G.add_edge(l1[p[0]], l2[p[1]])

# node_color = [node["color"] for node in G.nodes.values()]
# edge_color = [edge["color"] for edge in G.edges.values()]


nx.draw(G, pos, with_labels = False, node_color='gray', width=0.5)

plt.savefig("neuralnet.pdf")
