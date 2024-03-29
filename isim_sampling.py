from isim_comp import *
import math 

def medoid_sampling(fingerprints, n_ary = 'JT', percentage = 10):
    """
    This function samples a percentage of the objects with the lowest complementarity similarity, the medoids.
    
    Parameters:
    fingerprints: numpy array of fingerprints
    n_ary: type of similarity to index compute
    percentage: percentage of objects to sample
    
    Returns:
    indexes: indexes of the sampled objects
    """
    # Compute the complementarity similarity matrix
    comp_sim =  calculate_comp_sim(fingerprints, n_ary = n_ary)
    
    # Sort the complementarity similarities and get the indexes of the sorted array
    indexes = np.argsort(comp_sim)
    indexes = indexes[:int(len(indexes)*percentage/100)]
    
    return indexes

def outlier_sampling(fingerprints, n_ary = 'JT', percentage = 10):
    """
    This function samples a percentage of the objects with the highest complementarity similarity, the outliers.
    
    Parameters:
    fingerprints: numpy array of fingerprints
    n_ary: type of similarity to index compute
    percentage: percentage of objects to sample
    
    Returns:
    indexes: indexes of the sampled objects
    """
    # Compute the complementarity similarity matrix
    comp_sim =  calculate_comp_sim(fingerprints, n_ary = n_ary)
    
    # Sort the complementarity similarities and get the indexes of the sorted array
    indexes = np.argsort(comp_sim)
    indexes = indexes[-int(len(indexes)*percentage/100):]
    
    return indexes

def extremes_sampling(fingerprints, n_ary = 'JT', percentage = 10):
    """
    This function samples a percentage of the objects with the highest and lowest complementarity similarity, medoids and outliers.
    
    Parameters:
    fingerprints: numpy array of fingerprints
    n_ary: type of similarity to index compute
    percentage: percentage of objects to sample
    
    Returns:
    indexes: indexes of the sampled objects
    """
    # Define the percentage of extremes to sample
    percentage = percentage/2

    # Compute the complementarity similarity matrix
    comp_sim =  calculate_comp_sim(fingerprints, n_ary = n_ary)
    
    # Sort the complementarity similarities and get the indexes of the sorted array
    indexes = np.argsort(comp_sim)
    indexes = np.concatenate((indexes[:int(len(indexes)*percentage/100)], indexes[-int(len(indexes)*percentage/100):]))
    
    return indexes

def batched_sampling(fingerprints, n_ary = 'JT', percentage = 10, batches = None):
    """
    This function separates the objects in batches according to their complementarity similarity and samples a percentage of the objects
    in each batch to add up to the desired total percentage.
    
    Parameters:
    fingerprints: numpy array of fingerprints
    n_ary: type of similarity to index compute
    percentage: percentage of objects to sample
    batches: number of batches to separate the objects
    
    Returns:
    sampled_indexes: indexes of the sampled objects
    """
    # Define the number of objects
    n_objects = len(fingerprints)

    # Compute the complementarity similarity matrix
    comp_sim =  calculate_comp_sim(fingerprints, n_ary = n_ary)
    
    # Sort the complementarity similarities and get the indexes of the sorted array
    indexes = np.argsort(comp_sim)

    # Define the number of batches if not specified
    if not batches:
        batches = int(n_objects*percentage/100)

    # Split the data in batches
    batches = np.array_split(indexes, batches)    
        
    # Define the total number of objects to sample and the number of objects to sample in each batch
    n_sample = int(n_objects*percentage/100)

    # Check if the number of objects to sample is not less than the number of batches
    if n_sample < len(batches):
        raise ValueError("Warning: The number of objects to sample is too low for the number of batches, please specify a higher percentage, or a lower number of batches")

    # Sample the objects in each batch
    sampled_indexes = []
    i = 0
    while len(sampled_indexes) < n_sample:
        for b in batches:
            if len(b) > i:
                sampled_indexes.append(b[i])
                if len(sampled_indexes) >= n_sample:
                        break
            else:
                pass
        i += 1    
  
    return np.array(sampled_indexes)

def representative_sampling(fingerprints, n_ary = 'JT', percentage = 10, n_bins = 10, hard_cap = True):
    """
    Representative sampling according to comp_sim values.
    
    Divides the range of comp_sim values in nbins and then
    uniformly selects nsample molecules, consecutively
    taking one from each bin.
    
    Parameters:
    fingerprints: numpy array of fingerprints
    n_ary: type of similarity index to compute
    percentage: percentage of objects to sample
    
    Returns:
    sampled_indexes: indexes of the sampled objects
    """
    # Define the  number of objects and the number of objects to sample
    n_objects = len(fingerprints)
    n_sample = int(n_objects*percentage/100)

    # Check if the number of objects to sample is not less than the number of bins
    if n_sample < 1 or n_sample < n_bins:
        raise ValueError("Warning: The number of objects to sample is too low for the number of bins, please specify a higher percentage, or a lower number of bins")
    
    # Compute the complementarity similarity matrix
    comp_sim =  calculate_comp_sim(fingerprints, n_ary = n_ary)

    # Get the min and max comp_sim values
    min = np.min(comp_sim)
    max = np.max(comp_sim)

    # Divide the range of comp_sim values in n_bins
    D = max - min
    step = D/n_bins
    
    # Separate the objects in bins
    bins = []

    indices = np.array(list(range(n_objects)))
    for i in range(n_bins - 1):
        low = min + i * step
        up = min + (i + 1) * step
        ind = indices[(comp_sim >= low) * (comp_sim < up)]
        bin_comp_sim = comp_sim[ind]
        bins.append(ind[np.argsort(bin_comp_sim)])
    ind = indices[(comp_sim >= up) * (comp_sim <= max)]
    bin_comp_sim = comp_sim[ind]
    bins.append(ind[np.argsort(bin_comp_sim)])

    # Sample the objects from each bin
    order_sampled = []
    i = 0
    while len(order_sampled) < n_sample:
        for b in bins:
            if len(b) > i:
                order_sampled.append(b[i])
                if hard_cap:
                    if len(order_sampled) >= n_sample:
                        break
            else:
                pass
        i += 1

    return np.array(order_sampled)

