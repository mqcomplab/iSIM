import numpy as np
import random
import glob
import pickle
import math
from math import log
import time
import matplotlib.pyplot as plt

"""                         REAL-VALUED iSIM_MODULES
    ----------------------------------------------------------------------
    
    Miranda-Quintana Group, Department of Chemistry, University of Florida 
    
    ----------------------------------------------------------------------
    
    Please, cite the original papers on the n-ary indices:

    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00505-3
    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00504-4    
    """

def process_matrix(data):
    sq_data = data**2
    data_sum = np.sum(data, axis = 0)
    sq_data_sum = np.sum(sq_data, axis = 0)
    ij_array = 0.5 * (data_sum**2 - sq_data_sum)
    ij = np.sum(ij_array)
    return data_sum, sq_data_sum, ij

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

#total_max = 500
#size_max = 1000
#s = ''
#for _ in range(100):
#    fp_total = np.random.randint(10, total_max)
#    fp_size = np.random.randint(10, size_max)
    
#    total_fingerprints = np.random.rand(fp_total, fp_size)
    
    #total_fingerprints = np.array([[1,1,0,1,0],[1,0,1,1,0],[0,1,1,1,0],[0,0,0,1,1]])
    
    #isim = gen_sim_dict(total_fingerprints)['JT']
#    isim = isim_jt(total_fingerprints)
    
#    jts = []
#    for i, fp1 in enumerate(total_fingerprints):
#        for j, fp2 in enumerate(total_fingerprints):
#            if i == j:
#                pass
#            else:
                #jt = np.dot(fp1, fp2)/(np.dot(fp1, fp1)+np.dot(fp2, fp2)-np.dot(fp1, fp2))
                #jt2 = np.dot(fp1, fp2)/(np.dot(fp1, fp2) + np.dot(1-fp1, fp2) + np.dot(fp1, 1-fp2))
                #print(jt, jt2)
                #jts.append(pair_jt(fp1, fp2))
#    jts = np.array(jts)
#    p_jt = np.mean(jts)
    #print(pair_jt, isim)
#    s += '{:10.6}{:10.6}\n'.format(p_jt, isim)
#with open('test.txt', 'w') as outfile:
#    outfile.write(s[:-1])
