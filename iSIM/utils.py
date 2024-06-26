from iSIM.comp import calculate_isim
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors

def pairwise_average(fingerprints: np.ndarray, n_ary: str = 'RR'):
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

def binary_fps(smiles: list, fp_type: str = 'RDKIT', n_bits: int = 2048):
    """
    This function generates binary fingerprints for the dataset.
    
    Parameters:
    smiles: list of SMILES strings
    fp_type: type of fingerprint to generate ['RDKIT', 'ECFP4', 'ECFP6', or 'MACCS']
    n_bits: number of bits for the fingerprint
    
    Returns:
    fingerprints: numpy array of fingerprints
    """
    # Generate the fingerprints
    if fp_type == 'RDKIT':
       def generate_fp(mol, fp):
            DataStructs.cDataStructs.ConvertToNumpyArray(Chem.RDKFingerprint(mol), fp)
    elif fp_type == 'ECFP4':
        def generate_fp(mol, fp):
            DataStructs.cDataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits), fp)
    elif fp_type == 'ECFP6':
        def generate_fp(mol, fp):
            DataStructs.cDataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=n_bits), fp)
    elif fp_type == 'MACCS':
        def generate_fp(mol, fp):
            DataStructs.cDataStructs.ConvertToNumpyArray(Chem.MACCSkeys.GenMACCSKeys(mol), fp)
    else:
        print('Invalid fingerprint type: ', fp_type)
        exit(0)

    fingerprints = []
    for smi in smiles:
        # Generate the mol object
        try:
          mol = Chem.MolFromSmiles(smi)
        except:
          print('Invalid SMILES: ', smi)
          exit(0)

        # Generate the fingerprint and append to the list
        fingerprint = np.array([])
        generate_fp(mol, fingerprint)
        fingerprints.append(fingerprint)
    
    fingerprints = np.array(fingerprints)

    return fingerprints

def real_fps(smiles):
    """
    This function generates real number fingerprints for the dataset.
    
    Parameters:
    smiles: list of SMILES strings
    
    Returns:
    fingerprints: numpy array of fingerprints
    """
    fps = []
    for smi in smiles:
        # Generate the mol object
        try:
          mol = Chem.MolFromSmiles(smi)
        except:
          print('Invalid SMILES: ', smi)
          exit(0)

        # Generate the fingerprint and append to the list
        des = []
        for nm, fn in Descriptors._descList:
            try: 
                val = fn(mol)
            except:
                print('Error computing descriptor: ', nm)
                val = 'NaN'
            des.append(val)

        fps.append(des)
    
    # Drop columns with NaN values
    fps = np.array(fps)
    fps = fps[:, ~np.isnan(fps).any(axis = 0)]
    
    return fps

def minmax_norm(fps):
    """
    This function performs min-max normalization on the dataset.

    Parameters:
    fps: numpy array of fingerprints

    Returns:
    fps: normalized numpy array of fingerprints
    """

    # Turn the array into a DataFrame
    df = pd.DataFrame(fps)

    # Normalize the data
    df_numeric = df.select_dtypes(include = [np.number])
    columns = df_numeric.columns

    for column in columns:
        min_prop = np.min(df[column])
        max_prop = np.max(df[column])

        try:
            df[column] = [(x - min_prop) / (max_prop - min_prop) for x in df[column]]
        except ZeroDivisionError:
            df.drop(column, axis = 1)

    df = df.dropna(axis = 'columns')

    # Return the normalized data as a numpy array
    return df.to_numpy()

def rdkit_pairwise_sim(fingerprints):
    nfps = len(fingerprints)
    similarity = []

    for n in range(nfps - 1):
        sim = DataStructs.BulkTanimotoSimilarity(fingerprints[n], fingerprints[n+1:])
        similarity.extend([s for s in sim])

    return np.mean(similarity)
