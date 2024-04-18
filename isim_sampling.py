from isim_comp import *

def medoid_sampling(fingerprints = None, n_ary = 'JT', percentage = 10, comp_sim = None):
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
    if comp_sim is None:
        comp_sim =  calculate_comp_sim(fingerprints, n_ary = n_ary)
    else:
        comp_sim = comp_sim
    
    # Sort the complementarity similarities and get the indexes of the sorted array
    indexes = np.argsort(comp_sim)
    indexes = indexes[:int(len(indexes)*percentage/100)]
    
    return indexes

def outlier_sampling(fingerprints = None, n_ary = 'JT', percentage = 10, comp_sim = None):
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
    if comp_sim is None:
        comp_sim =  calculate_comp_sim(fingerprints, n_ary = n_ary)
    else:
        comp_sim = comp_sim

    # Sort the complementarity similarities and get the indexes of the sorted array
    indexes = np.argsort(comp_sim)
    indexes = indexes[-int(len(indexes)*percentage/100):]
    
    return indexes

def extremes_sampling(fingerprints = None, n_ary = 'JT', percentage = 10, comp_sim = None):
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

    # Compute the complementarity similarity matrix, changes to see if comp_sim is provided
    if comp_sim is None:
        comp_sim =  calculate_comp_sim(fingerprints, n_ary = n_ary)
    else:
        comp_sim = comp_sim
    
    # Sort the complementarity similarities and get the indexes of the sorted array
    indexes = np.argsort(comp_sim)
    indexes = np.concatenate((indexes[:int(len(indexes)*percentage/100)], indexes[-int(len(indexes)*percentage/100):]))
    
    return indexes

def stratified_sampling(fingerprints = None, n_ary = 'JT', percentage = 10, strata = None, comp_sim = None):
    """
    This function separates the objects in strata according to their complementarity similarity and it samples a percentage of the objects
    in each batch to add up to the desired total percentage. If objects to sample in each stratum are not equal lowest complementary
    similarity strata are sampled first.
    
    Parameters:
    fingerprints: np.ndarray, numpy array of fingerprints
    n_ary: str, type of similarity to index compute {'JT', 'SM', 'RR'}
    percentage: int or float, percentage of objects to sample
    strata: int, number of strata to separate the objects in
    
    Returns:
    sampled_indexes: indexes of the sampled objects
    """
    # Compute the complementarity similarity matrix
    if comp_sim is None:
        comp_sim =  calculate_comp_sim(fingerprints, n_ary = n_ary)
        n_objects = len(fingerprints)
    else:
        comp_sim = comp_sim
        n_objects = len(comp_sim)
    
    # Sort the complementarity similarities and get the indexes of the sorted array
    indexes = np.argsort(comp_sim)

    # Define the number of batches if not specified
    if not strata:
        strata = int(n_objects*percentage/100)

    # Split the data in batches
    strata = np.array_split(indexes, strata)    
        
    # Define the total number of objects to sample and the number of objects to sample in each batch
    n_sample = int(n_objects*percentage/100)

    # Check if the number of objects to sample is not less than the number of batches
    if n_sample < len(strata):
        raise ValueError("Warning: The number of objects to sample is too low for the number of batches, please specify a higher percentage, or a lower number of batches")

    # Sample the objects in each batch
    sampled_indexes = []
    i = 0
    while len(sampled_indexes) < n_sample:
        for stratum in strata:
            if len(stratum) > i:
                sampled_indexes.append(stratum[i])
                if len(sampled_indexes) >= n_sample:
                        break
            else:
                pass
        i += 1    
  
    return np.array(sampled_indexes)

def quota_sampling(fingerprints = None, n_ary = 'JT', percentage = 10, n_bins = 10, hard_cap = True, comp_sim = None):
    """
    Quota sampling according to comp_sim values.
    
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
    
    # Compute the complementarity similarity matrix
    if comp_sim is None:
        comp_sim =  calculate_comp_sim(fingerprints, n_ary = n_ary)
        n_objects = len(fingerprints)
    else:
        comp_sim = comp_sim
        n_objects = len(comp_sim)
    
    # Define the number of objetcs to sample
    n_sample = int(n_objects*percentage/100)
    
    # Check if the number of objects to sample is not less than the number of bins
    if n_sample < 1 or n_sample < n_bins:
        raise ValueError("Warning: The number of objects to sample is too low for the number of bins, please specify a higher percentage, or a lower number of bins")
    
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