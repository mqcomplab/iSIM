from isim_comp import gen_sim_dict
import numpy as np


""" Proof-of-principle Hierarchical Agglomerative Clustering using iSIM
    ----------------------------------------------------------------------
    
    Miranda-Quintana Group, Department of Chemistry, University of Florida 
    
    ----------------------------------------------------------------------
    
    Please, cite the original paper on iSIM:

    https://pubs.rsc.org/en/content/articlehtml/2024/dd/d4dd00041b

    """

def pre_process(fingerprints):
    """
    Function to pre-process the fingerprints before clustering.
    
    Parameters:
    -----------
    fingerprints: numpy.ndarray or list
        Array of fingerprints to cluster.
    
    Returns:
    --------
    processed_fingerprints: numpy.ndarray
        Array of pre-processed fingerprints.
    """
    processed_fingerprints = []
    for i, fp in enumerate(fingerprints):
        l = np.append(fp, 1)
        processed_fingerprints.append([np.array(l), [i]])

    return processed_fingerprints

def max_indices(fingerprints, n_ary = 'RR'):
    """
    Function to find the indices of the two most similar fingerprints in a list of fingerprints.

    Parameters:
    -----------
    fingerprints: numpy.ndarray or list
        Array of fingerprints to cluster.
    n_ary: str
        Type of iSIM to use. Default is 'RR'. Other options are 'JT' and 'SM'.

    Returns:
    --------
    max1: int
        Index of the first most similar fingerprint.
    max2: int
        Index of the second most similar fingerprint.
    """
    max_sim = -3.08
    max1, max2 = len(fingerprints), len(fingerprints)
    for i, d1 in enumerate(fingerprints):
        for j, d2 in enumerate(fingerprints):
            if i == j:
                pass
            else:
                fp1 = np.array(d1[0][:-1])
                fp2 = np.array(d2[0][:-1])
                n = d1[0][-1] + d2[0][-1]
                s = gen_sim_dict(data = fp1 + fp2, n_objects = n)[n_ary]
                if s > max_sim:
                    max_sim = s
                    max1 = i
                    max2 = j
    return max1, max2

def update_data(data, max1, max2):
    """
    Function to update the list of fingerprints after combining the two most similar fingerprints.

    Parameters:
    -----------
    data: list
        List of fingerprints to cluster.
    max1: int
        Index of the first most similar fingerprint.
    max2: int
        Index of the second most similar fingerprint.
    
    Returns:
    --------
    new_data: list
        List of fingerprints after combining the two most similar fingerprints.
    """
    new_cluster = data[max1][-1] + data[max2][-1]
    new_condensed = data[max1][0] + data[max2][0]
    new_data = []
    for i, d in enumerate(data):
        if i == max1:
            pass
        elif i == max2:
            pass
        else:
            new_data.append(d)
    new_data.append([new_condensed, new_cluster])
    return new_data

def cluster_tree(processed_fingerprints, n_ary = 'RR'):
    tree = []
    while len(processed_fingerprints) > 1:
        tree.append([])
        for d in processed_fingerprints:
            tree[-1].append(d[-1])
        max1, max2 = max_indices(processed_fingerprints, n_ary = n_ary)
        processed_fingerprints = update_data(processed_fingerprints, max1, max2)
    return tree

def gen_z(tree):
    # Combine the last two elements of the tree
    tree.append([tree[-1][-1] + tree[-1][-2]])

    # Get the numer of original data objects
    n = len(tree[0])

    # Initialize the Z matrix
    Z = np.zeros((n - 1, 4))

    # Fill the third and fourth elements of the Z matrix
    for i in range(0, n-1):
        Z[i, 2] = i + 1
        Z[i, 3] = len(tree[i + 1][-1])

    # Get the clusters in order for indexing
    clusters = []
    for item in tree[0]:
        clusters.append(item)

    for i in range(1, len(tree) - 1):
        clusters.append(tree[i][-1])
    
    # Fill the first and second elements of the Z matrix
    for i in range(0, n - 1):
        # Get elements in the cluster
        c_clusters = []
        for element in tree[i]:
            if element not in tree[i+1]:
                c_clusters.append(element)
    
        # Combine elements in c_clusters in one list
        comb_clusters = [item for sublist in c_clusters for item in sublist]

        #if comb_clusters == tree[i+1][-1] and len(c_clusters) == 2: print('OK')
        Z[i][0] = clusters.index(c_clusters[0])
        Z[i][1] = clusters.index(c_clusters[1])

    return Z

def hierarchical_clustering(fingerprints, n_ary = 'RR'):
    processed_fingerprints = pre_process(fingerprints)
    tree = cluster_tree(processed_fingerprints, n_ary = n_ary)
    Z = gen_z(tree)
    return tree, Z

