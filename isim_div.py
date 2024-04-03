from isim_comp import *
import random

def get_new_index_n(total_data, selected_condensed, n, select_from_n, n_ary = 'RR'):
    """Select a diverse object using the ECS_MeDiv algorithm"""
    n_total = n + 1

    # min value that is guaranteed to be higher than all the comparisons
    min_value = 3.08
    
    # placeholder index
    index = len(total_data) + 1
    
    # for all indices that have not been selected
    for i in select_from_n:
        # column sum
        c_total = selected_condensed + total_data[i]
        # calculating similarity
        sim_index = gen_sim_dict(c_total, n_total)[n_ary]
        # if the sim of the set is less than the similarity of the previous diverse set, update min_value and index
        if sim_index < min_value:
            index = i
            min_value = sim_index
    return index

def get_new_index_sqrt(total_data, selected_condensed, n, select_from_n, k = 2, n_ary = 'RR'):
    n_total = n + 1

    # min value that is guaranteed to be higher than all the comparisons
    min_value = 3.08
    
    # placeholder index
    index = len(total_data) + 1
    
    # for all indices that have not been selected
    for i in select_from_n:
        # column sum
        c_total = selected_condensed + total_data[i]
        # calculating similarity
        sim_index = gen_sim_dict(c_total, n_objects = n_total, k = k)[n_ary]
        # if the sim of the set is less than the similarity of the previous diverse set, update min_value and index
        if sim_index < min_value:
            index = i
            min_value = sim_index
    return index

def get_new_index_reverse(total_data, selected_condensed, n, select_from_n, n_ary = 'RR'):
    n_total =  n - 1

    #min value that is guaranteed to be lower than all the comparisons
    min_value = 3.08

    #placeholder index 
    index = len(total_data) + 1

    #for all the molecules that are already selected
    for i in select_from_n:
        #new column sum
        c_total = selected_condensed - total_data[i]
        #calculating isim
        sim_index = gen_sim_dict(c_total, n_objects = n_total)[n_ary]
        # if the sim of the set when taking that molecule out is lower than the min value, store new index and value
        if sim_index < min_value:
            index = i
            min_value = sim_index

    return index


def get_new_indices_b_max(total_data, selected_b, select_from_b, n_ary):
    all_comps = []
    for i in select_from_b:
        comps = []
        for j in selected_b:
            new_indices = [i, j]
            new_fingerprints = total_data[new_indices]
            sim_index = gen_sim_dict(new_fingerprints)[n_ary]
            comps.append(sim_index)
        all_comps.append(comps)
    sim_values = [max(comps) for comps in all_comps]
    min_sim = min(sim_values)
    min_list = [j for j, v in enumerate(sim_values) if v == min_sim]
    return select_from_b[min_list[0]]

def diversity(data, percentage: int, start = 'medoid', n_ary = 'RR', method = 'isim'):
    """ diversity: function to select from a dataset the most diverse molecules
    -----------------------------------------------------------------------

    Arguments
    ---------
    data: np.array
        Array of arrays containing the binary string objects 
     
    percentaje: int
        Percentage of the provided data that wants to be sampled

    start: str or list
        srt: key on what is used to start the selection 
        {'medoid', 'random', 'outlier'}  
         
        list: contains the indexes of the molecules you want to start the selection

    n_ary: str
        Key with the abbreviation of the similarity index to perform the selection 
    """
 
    # total number of objects
    n_total = len(data)

    # indices of all the objects
    total_indices = np.array(range(n_total))

    if start =='medoid':
        seed = calculate_medoid(data,  n_ary = n_ary)
        selected_n = [seed]
    elif start == 'random':
        seed = random.randint(0, n_total - 1)
        selected_n = [seed]
    elif start == 'outlier':
        seed = calculate_outlier(data, n_ary = n_ary)
        selected_n = [seed]
    elif isinstance(start, list):
        selected_n = start
    else:
        raise ValueError('Select a correct starting point: medoid, random or outlier')
  
 
    # Number of initial objects
    n = len(selected_n)

    # Number of objects be selected
    n_max = int(n_total * percentage / 100)

	# Condensation of selected initial selection 
    selected_condensed = np.sum([data[i] for i in selected_n], axis = 0)  

    while len(selected_n) < n_max:
        # indices from which to select the new fingerprints
        select_from_n = np.delete(total_indices, selected_n)

        if method == 'isim':
            # new index selected
            new_index_n = get_new_index_n(data, selected_condensed, n, select_from_n, n_ary = n_ary)
        elif method == 'bmax':
            new_index_n = get_new_indices_b_max(data, selected_n, select_from_n, n_ary = n_ary)
        elif isinstance(method, int):
            new_index_n = get_new_index_sqrt(data, selected_condensed, n, select_from_n, k = method, n_ary = n_ary)

        # updating column sum vector
        selected_condensed += data[new_index_n]

        # updating selected indices
        selected_n.append(new_index_n)
        n = len(selected_n)

    return selected_n

def reverse_diversity(data, percentage: int, n_ary = 'RR'):
    """ diversity: function to select from a dataset the most diverse molecules
    -----------------------------------------------------------------------

    Arguments
    ---------
    data: np.array
        Array of arrays containing the binary string objects 
     
    percentaje: int
        Percentage of the provided data that wants to be sampled

    n_ary: str
        Key with the abbreviation of the similarity index to perform the selection 
    """
 
    # total number of objects
    n_total = len(data)

    # indices of all the objects
    total_indices = np.array(range(n_total))

    # Number of initial objects
    n = n_total

    # Number of objects be deselected
    n_max = int(n_total * percentage / 100)

    # Deselected objeccts 
    deselected_n = []

    # Condensation of the total fingerprints 
    selected_condensed = np.sum(data, axis = 0)   

    # Select from list 
    select_from_n = total_indices 

    while len(select_from_n) > n_max:
        # indices from which to select the new fingerprints
        select_from_n = np.delete(total_indices, deselected_n)
        #print(selected_condensed)        

        # new index selected
        new_index_n = get_new_index_reverse(data, selected_condensed, n, select_from_n, n_ary = n_ary)

        # updating column sum vector
        selected_condensed = selected_condensed - data[new_index_n]

        # updating selected indices
        deselected_n.append(new_index_n)

        # Update the number of selected 
        n = n_total - len(deselected_n)

    return select_from_n


