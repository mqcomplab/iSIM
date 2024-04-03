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
