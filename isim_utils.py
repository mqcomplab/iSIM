from isim_comp import calculate_isim
from isim_real import *
import numpy as np

def pairwise_average(fingerprints, n_ary = 'RR'):
    """
    This function computes the pairwise average similarity between all objects in the dataset.
    
    Parameters:
    fingerprints: numpy array of fingerprints
    n_ary: type of similarity to index compute
    
    Returns:
    average: average similarity between all objects
    """
    # Compute the pairwise similarities
    pairwise_sims = []
    for i in range(len(fingerprints)):
        for j in range(len(fingerprints)):
            if i != j:
                pairwise_sims.append(calculate_isim(np.array([fingerprints[i], fingerprints[j]]), n_ary = n_ary))

    # Compute the average similarity
    average = np.mean(pairwise_sims)
    
    return average

def pairwise_average_real(fingerprints, n_ary = 'RR'):
    """
    This function computes the pairwise average similarity between all objects in the dataset for real numbers.
    
    Parameters:
    fingerprints: numpy array of fingerprints
    n_ary: type of similarity to index compute
    
    Returns:
    average: average similarity between all objects
    """
    if n_ary == 'RR': pair_func = pair_rr
    elif n_ary == 'SM': pair_func = pair_sm
    elif n_ary == 'JT': pair_func = pair_jt

    # Compute the pairwise similarities
    pairwise_sims = []
    for i in range(len(fingerprints)):
        for j in range(len(fingerprints)):
            if i != j:
                pairwise_sims.append(pair_func(fingerprints[i], fingerprints[j]))

    # Compute the average similarity
    average = np.mean(pairwise_sims)
    
    return average