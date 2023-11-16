import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from isim_comp import *

""" Proof-of-principle Hierarchical Agglomerative Clustering using iSIM
    ----------------------------------------------------------------------
    
    Miranda-Quintana Group, Department of Chemistry, University of Florida 
    
    ----------------------------------------------------------------------
    
    Please, cite the original paper on iSIM:

    """

# Create a cluster class
class Cluster:
    def __init__(self, index, c_sum, size, isim):
        self.index = index
        self.c_sum = c_sum
        self.isim = isim
        self.size = size

# Find the most similar clusters
def combine_clusters(clusters):
    max_isim = -3.08
    for i, cluster in enumerate(clusters):
        for j, cluster in enumerate(clusters):
            if i == j:
                pass
            else:
                c_sum = np.sum([clusters[i].c_sum, clusters[j].c_sum], axis = 0)
                size = clusters[i].size + clusters[j].size
                isim = gen_sim_dict(c_sum, size)['RR']

                if isim > max_isim:
                    max_isim = isim
                    new_cluster = Cluster(index = clusters[i].index + clusters[j].index, c_sum = c_sum, size = size, isim = isim)
                    del_index = [i, j]

    clusters = [cluster for i, cluster in enumerate(clusters) if i not in del_index]
    clusters.append(new_cluster)

    return clusters

# Combine clusters until there is only one cluster left
def hierarchical_clustering(fingerprints):
    clusters = []
    for i, fp in enumerate(fingerprints):
        clusters.append(Cluster(index = [i], c_sum = fp, size = 1, isim = 0))

    tree = []
    tree.append([cluster.index for cluster in clusters])

    while len(clusters) > 1:
        clusters = combine_clusters(clusters)
        print(len(clusters))
        tree.append([cluster.index for cluster in clusters])
    return tree