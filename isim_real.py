import numpy as np

"""                         REAL-VALUED iSIM_MODULES
    ----------------------------------------------------------------------
    
    Miranda-Quintana Group, Department of Chemistry, University of Florida 
    
    ----------------------------------------------------------------------
    
    Please, cite the original paper:   
    """

"""                         REAL-VALUED PAIRWISE SIMILARITY FUNCTIONS                                 """

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

def pair_jt(fp1, fp2):
    return np.dot(fp1, fp2)/(np.dot(fp1, fp1) + np.dot(fp2, fp2) - np.dot(fp1, fp2))

def pair_rr(fp1, fp2, m = None):
    if not m:
        m = len(fp1)
    return np.dot(fp1, fp2)/m

def pair_sm(fp1, fp2, m = None):
    if not m:
        m = len(fp1)
    return (np.dot(fp1, fp2) + np.dot(1 - fp1, 1 - fp2))/m


"""                         REAL-VALUE iSIM FUNCTIONS                                 """

def process_matrix(fingerprints: np.ndarray):
    """This function gives the necessary values for iSIM calculations for real-number fingerprints.
    Parameters:
    ----------------
    fingerprints: np.ndarray
        numpy array of fingerprints with real-valued numbers
    
    Returns:
    ----------------
    data_sum: np.ndarray
        columnwise sum of the fingerprints
    sq_data_sum: np.ndarray
        columnwise sum of the squared fingerprints
    ij: float
        sum of the inner products of the fingerprints"""

    # Calculate the Hadamard product of the fingerprints
    sq_data = fingerprints**2

    # Calculate the columnwise sum of the fingerprints and the squared fingerprints
    data_sum = np.sum(fingerprints, axis = 0)
    sq_data_sum = np.sum(sq_data, axis = 0)

    # Calculate linearly inner product of the fingerprints
    ij_array = 0.5 * (data_sum**2 - sq_data_sum)
    ij = np.sum(ij_array)

    return data_sum, sq_data_sum, ij

def calculate_isim_real(fingerprints: np.ndarray, n_ary: str = 'JT'):
    """This function calculates the iSIM index for real-valued fingerprints.
    Parameters:
    ----------------
    fingerprints: np.ndarray
        numpy array of fingerprints with real-valued numbers
    n_ary: str
        type of similarity index to calculate [JT, RR, SM]
    
    Returns:
    ----------------
    iSIM: float
        iSIM (average similarity)"""
    
    if n_ary == 'RR':
        return isim_rr(fingerprints)
    elif n_ary == 'SM':
        return isim_sm(fingerprints)
    elif n_ary == 'JT':
        return isim_jt(fingerprints)

def isim_jt(data):
    n_objects = len(data)
    data_sum, sq_data_sum, ij = process_matrix(data)
    inners = (n_objects - 1) * np.sum(sq_data_sum)
    return ij/(inners - ij)

def isim_rr(data):
    n_objects = len(data)
    m = len(data[0])
    data_sum, sq_data_sum, ij = process_matrix(data)
    return 2 * ij/(m * n_objects * (n_objects - 1))

def isim_sm(data):
    n_objects = len(data)
    m = len(data[0])
    data_sum, sq_data_sum, ij = process_matrix(data)
    flip_data = 1 - data
    flip_sum, sq_flip_sum, flip_ij = process_matrix(flip_data)
    return 2 * (ij + flip_ij)/(m * n_objects * (n_objects - 1))

"""                         REAL-VALUED COMPLEMENTARY SIMILARITY FUNCTIONS                                 """

def calculate_comp_sim_real(fingerprints: np.ndarray, n_ary: str = 'JT'):
    """This function calculates the composite similarity for real-valued fingerprints.
    Parameters:
    ----------------
    fingerprints: np.ndarray
        numpy array of fingerprints with real-valued numbers
    n_ary: str
        type of similarity index to calculate [JT, RR, SM]
    
    Returns:
    ----------------
    comp_sim: np.ndarray
        complementary similarity vector"""

    if n_ary == 'JT':
        comp_sim = comp_sim_jt(fingerprints)
    elif n_ary == 'RR':
        comp_sim = comp_sim_rr(fingerprints)
    elif n_ary == 'SM':
        comp_sim = comp_sim_sm(fingerprints)
    
    return comp_sim

def comp_sim_rr(fingerprints: np.ndarray):
    n = len(fingerprints) - 1
    had_matrix = fingerprints**2
    c_sum = np.sum(fingerprints, axis=0)
    sum_sq = np.sum(had_matrix, axis=0)
    comp_sq_sum = np.sum((c_sum - fingerprints)**2, axis=1)
    comp_sum_sq = np.sum(sum_sq - had_matrix, axis=1)

    comp_sim = np.array((comp_sq_sum - comp_sum_sq)/(len(fingerprints[0]) * n * (n - 1)))

    return comp_sim

def comp_sim_sm(fingerprints: np.ndarray):
    n = len(fingerprints) - 1
    had_matrix = fingerprints**2
    flip_matrix = 1 - fingerprints
    had_flip_matrix = flip_matrix**2
    c_sum = np.sum(fingerprints, axis=0)
    had_sum = np.sum(had_matrix, axis=0)
    flip_sum = np.sum(flip_matrix, axis=0)
    had_flip_sum = np.sum(had_flip_matrix, axis=0)

    comp_sq_sum = np.sum((c_sum - fingerprints)**2, axis=1)
    comp_sum_sq = np.sum(had_sum - had_matrix, axis=1)
    flip_sq_sum = np.sum((flip_sum - flip_matrix)**2, axis=1)
    flip_sum_sq = np.sum(had_flip_sum - had_flip_matrix, axis=1)

    comp_sim = np.array((comp_sq_sum - comp_sum_sq + flip_sq_sum - flip_sum_sq)/(len(fingerprints[0]) * n * (n - 1)))
    
    return comp_sim

def comp_sim_jt(fingerprints: np.ndarray):
    n = len(fingerprints) - 1
    c_sum = np.sum(fingerprints, axis=0)
    sum_sq = np.sum(fingerprints**2, axis=0)
    comp_sq_sum = np.sum((c_sum - fingerprints)**2, axis=1)
    comp_sum_sq = np.sum(sum_sq - fingerprints**2, axis=1)

    comp_sim = np.array(0.5*(comp_sq_sum - comp_sum_sq)/((n - 1)*comp_sum_sq-0.5*(comp_sq_sum - comp_sum_sq)))

    return comp_sim

def calculate_medoid_real(data, n_ary = 'JT'):
    return np.argmin(calculate_comp_sim_real(data, n_ary = n_ary))

def calculate_outlier_real(data, n_ary = 'JT'):
    return np.argmax(calculate_comp_sim_real(data, n_ary = n_ary))