# iSIM
Module containing scripts to perform multiple comparisons simultaneously and getting the exact same value as the average pairwise comparisons of molecules represented by binary fingerprints or real number descriptors.

Pre-print
https://chemrxiv.org/engage/chemrxiv/article-details/65567ce56e0ec7777f17b138

The quantification of molecular similarity has been present since the beginning of cheminformatics. Although several similarity indices and molecular representations have been reported, all of them ultimately reduce to the calculation of molecular similarities of only two objects at a time. Hence, to get the average similarity of a set of molecules, all the pairwise comparisons need to be computed, which demands a quadratic scaling in the number of computational resources. Here we propose an exact alternative to this problem: iSIM (Instant Similarity). iSIM performs comparisons of multiple molecules at the same time and yields the same value as the average pairwise comparisons of molecules represented with binary fingerprints and real-value descriptors. In this work, we introduce the mathematical framework and several applications of iSIM in chemical sampling, visualization, diversity selection, and clustering.


Curated CHEMBL datasets were retrieved from:
https://github.com/molML/MoleculeACE/blob/main/MoleculeACE/Data/benchmark_data/old/
"""Exposing the Limitations of Molecular Machine Learning with Activity Cliffs. Derek van Tilborg, Alisa Alenicheva, and Francesca Grisoni. Journal of Chemical Information and Modeling, 2022, 62 (23), 5938-5951. DOI: 10.1021/acs.jcim.2c01073"""
