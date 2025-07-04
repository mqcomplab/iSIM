import numpy as np

"""                             iSIM_MODULES
    ----------------------------------------------------------------------
    
    Miranda-Quintana Group, Department of Chemistry, University of Florida 
    
    ----------------------------------------------------------------------
    
    Please, cite the original paper on iSIM:

    https://pubs.rsc.org/en/content/articlehtml/2024/dd/d4dd00041b

    """

def calculate_counters(data, n_objects = None, k = 1):
    """Calculate 1-similarity, 0-similarity, and dissimilarity counters

    Arguments
    ---------
    data : np.ndarray
        Array of arrays, each sub-array contains the binary object 
        OR Array with the columnwise sum, if so specify n_objects

    n_objects : int
        Number of objects, only necessary if the column wize sum is the input data.
    
    k : int
        Integer indicating the 1/k power used to approximate the average of the
        similarity values elevated to 1/k.

    Returns
    -------
    counters : dict
        Dictionary with the weighted and non-weighted counters.

    """
    
    # Check if the data is a np.ndarray of a list
    if not isinstance(data, np.ndarray):
        raise TypeError("Warning: Input data is not a np.ndarray, to secure the right results please input the right data type")
    
    if data.ndim == 1:
        c_total = data
        if not n_objects:
            raise ValueError("Input data is the columnwise sum, please specify number of objects")
    else:
        c_total = np.sum(data, axis = 0)
        if not n_objects:
            n_objects = len(data)      
        elif n_objects and n_objects != len(data):
            print("Warning, specified number of objects is different from the number of objects in data")
            n_objects = len(data)
            print("Doing calculations with", n_objects, "objects.")

    # Calculate a, d, b + c
    a_array = c_total * (c_total - 1) / 2
    off_coincidences = n_objects - c_total
    d_array = off_coincidences * (off_coincidences - 1) / 2
    dis_array = off_coincidences * c_total

    a = np.sum(np.power(a_array, 1/k))
    d = np.sum(np.power(d_array, 1/k))
    total_dis = np.sum(np.power(dis_array, 1/k))
            
    total_sim = a + d
    p = total_sim + total_dis
    
    counters = {"a": a, "d": d, "total_sim": total_sim,
                "total_dis": total_dis, "p": p}
    return counters

def calculate_isim(data, n_objects = None, n_ary = 'RR'):
    """Calculate the iSIM index for RR, JT, or SM

    Arguments
    ---------
    data : np.ndarray
        Array of arrays, each sub-array contains the binary object 
        OR Array with the columnwise sum, if so specify n_objects
    
    n_objects : int
        Number of objects, only necessary if the column wize sum is the input data.

    n_ary : str
        String with the initials of the desired similarity index to calculate the iSIM from. 
        Only RR, JT, or SM are available. For other indexes use gen_sim_dict.

    Returns
    -------
    isim : float
        iSIM index for the specified similarity index.
    """

    # Check if the data is a np.ndarray of a list
    if not isinstance(data, np.ndarray):
        raise TypeError("Warning: Input data is not a np.ndarray, to secure the right results please input the right data type")
    
    if data.ndim == 1:
        c_total = data
        if not n_objects:
            raise ValueError("Input data is the columnwise sum, please specify number of objects")
    else:
        c_total = np.sum(data, axis = 0)
        if not n_objects:
            n_objects = len(data)      
        elif n_objects and n_objects != len(data):
            print("Warning, specified number of objects is different from the number of objects in data")
            n_objects = len(data)
            print("Doing calculations with", n_objects, "objects.")

    # Calculate only necessary counters for the desired index 

    if n_ary == 'RR':
        a = np.sum(c_total * (c_total - 1) / 2)
        p = n_objects * (n_objects - 1) * len(c_total) / 2

        return a/p
    
    elif n_ary == 'JT':
        a = np.sum(c_total * (c_total - 1) / 2)
        off_coincidences = n_objects - c_total
        total_dis = np.sum(off_coincidences * c_total)

        return a/(a + total_dis)
    
    elif n_ary == 'SM':
        a = np.sum(c_total * (c_total - 1) / 2)
        off_coincidences = n_objects - c_total
        d = np.sum(off_coincidences * (off_coincidences - 1) / 2)
        p = n_objects * (n_objects - 1) * len(c_total) / 2

        return (a + d)/p

def gen_sim_dict(data, n_objects = None, k = 1):
    """Calculate a dictionary containing all the available similarity indexes

    Arguments
    ---------
    See calculate counters.

    Returns
    -------
    sim_dict : dict
        Dictionary with the weighted and non-weighted similarity indexes."""
 
    # Indices
    # AC: Austin-Colwell, BUB: Baroni-Urbani-Buser, CTn: Consoni-Todschini n
    # Fai: Faith, Gle: Gleason, Ja: Jaccard, Ja0: Jaccard 0-variant
    # JT: Jaccard-Tanimoto, RT: Rogers-Tanimoto, RR: Russel-Rao
    # SM: Sokal-Michener, SSn: Sokal-Sneath n
 
    # Calculate the similarity and dissimilarity counters
    counters = calculate_counters(data = data, n_objects = n_objects, k = k)

    ac = (2/np.pi) * np.arcsin(np.sqrt(counters['total_sim']/
                                       counters['p']))
    bub = ((counters['a'] * counters['d'])**0.5 + counters['a'])/\
          ((counters['a'] * counters['d'])**0.5 + counters['a'] + counters['total_dis'])
    fai = (counters['a'] + 0.5 * counters['d'])/\
          (counters['p'])
    gle = (2 * counters['a'])/\
          (2 * counters['a'] + counters['total_dis'])
    ja = (3 * counters['a'])/\
         (3 * counters['a'] + counters['total_dis'])
    jt = (counters['a'])/\
         (counters['a'] + counters['total_dis'])
    rt = (counters['total_sim'])/\
         (counters['p'] + counters['total_dis'])
    rr = (counters['a'])/\
         (counters['p'])
    sm = (counters['total_sim'])/\
         (counters['p'])
    ss1 = (counters['a'])/\
          (counters['a'] + 2 * counters['total_dis'])
    ss2 = (2 * counters['total_sim'])/\
          (counters['p'] + counters['total_sim'])

    # Dictionary with all the results
    Indices = {'AC': ac, 'BUB':bub, 'Fai':fai, 'Gle':gle, 'Ja':ja,
               'JT':jt, 'RT':rt, 'RR':rr, 'SM':sm, 'SS1':ss1, 'SS2':ss2}
    #Indices = {'Fai':fai, 'Gle':gle, 'Ja':ja,
    #           'JT':jt, 'RT':rt, 'RR':rr, 'SM':sm, 'SS1':ss1, 'SS2':ss2}
    return Indices

def calculate_medoid(data, n_ary = 'RR'):
    return np.argmin(calculate_comp_sim(data, n_ary = n_ary))

def calculate_outlier(data, n_ary = 'RR'):
    return np.argmax(calculate_comp_sim(data, n_ary = n_ary))


def calculate_comp_sim(data, n_ary = 'RR'):
    """Calculate the complementary similarity for RR, JT, or SM

    Arguments
    ---------
    data : np.ndarray
        Array of arrays, each sub-array contains the binary object 
        
    n_objects : int
        Number of objects, only necessary if the column wize sum is the input data.

    n_ary : str
        String with the initials of the desired similarity index to calculate the iSIM from. 
        Only RR, JT, or SM are available. For other indexes use gen_sim_dict.

    Returns
    -------
    comp_sims : nd.array
        1D array with the complementary similarities of all the molecules in the set.
    """
    
    colsum = np.sum(data, axis = 0)
    n_objects = len(data) - 1   
    comp_sims = [calculate_isim(colsum - data[i], n_objects = n_objects, n_ary = n_ary) for i in range(len(data))]
    
    return comp_sims