# iSIM
Module containing scripts to perform multiple comparisons simultaneously and getting the exact same value as the average pairwise comparisons of molecules represented by binary fingerprints or real number descriptors.

Paper
https://pubs.rsc.org/en/content/articlehtml/2024/dd/d4dd00041b

The quantification of molecular similarity has been present since the beginning of cheminformatics. Although several similarity indices and molecular representations have been reported, all of them ultimately reduce to the calculation of molecular similarities of only two objects at a time. Hence, to get the average similarity of a set of molecules, all the pairwise comparisons need to be computed, which demands a quadratic scaling in the number of computational resources. Here we propose an exact alternative to this problem: iSIM (Instant Similarity). iSIM performs comparisons of multiple molecules at the same time and yields the same value as the average pairwise comparisons of molecules represented with binary fingerprints and real-value descriptors. In this work, we introduce the mathematical framework and several applications of iSIM in chemical sampling, visualization, diversity selection, and clustering.


Curated CHEMBL datasets were retrieved from:
https://github.com/molML/MoleculeACE/blob/main/MoleculeACE/Data/benchmark_data/old/
"""Exposing the Limitations of Molecular Machine Learning with Activity Cliffs. Derek van Tilborg, Alisa Alenicheva, and Francesca Grisoni. Journal of Chemical Information and Modeling, 2022, 62 (23), 5938-5951. DOI: 10.1021/acs.jcim.2c01073"""

### To install: 
`git clone git@github.com:mqcomplab/iSIM.git`

`pip install -e .`
### Prerequisites
iSIM current supports Python 3.10. 

For iSIM functions:

numpy (1.24.2)

For plots and analysis (notebooks):

matplotlib (3.7.3)

pandas (2.0.2)

random 

rdkit (2023.03.1)

scipy (1.10.1)

seaborn (0.12.2)

sklearn (1.3.1)

### Funding
Research contained in this package was supported by the National Institute of General Medical Sciences of the National Institutes  of  Health  under  award  number  R35GM150620. 

