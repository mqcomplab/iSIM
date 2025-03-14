# iSIM
Module containing scripts to perform multiple comparisons simultaneously and getting the exact same value as the average pairwise comparisons of molecules represented by binary fingerprints or real number descriptors.

Paper
https://pubs.rsc.org/en/content/articlehtml/2024/dd/d4dd00041b

The quantification of molecular similarity has been present since the beginning of cheminformatics. Although several similarity indices and molecular representations have been reported, all of them ultimately reduce to the calculation of molecular similarities of only two objects at a time. Hence, to get the average similarity of a set of molecules, all the pairwise comparisons need to be computed, which demands a quadratic scaling in the number of computational resources. Here we propose an exact alternative to this problem: iSIM (Instant Similarity). iSIM performs comparisons of multiple molecules at the same time and yields the same value as the average pairwise comparisons of molecules represented with binary fingerprints and real-value descriptors. In this work, we introduce the mathematical framework and several applications of iSIM in chemical sampling, visualization, diversity selection, and clustering.


Curated CHEMBL datasets were retrieved from:
https://github.com/molML/MoleculeACE/blob/main/MoleculeACE/Data/benchmark_data/old/
"""Exposing the Limitations of Molecular Machine Learning with Activity Cliffs. Derek van Tilborg, Alisa Alenicheva, and Francesca Grisoni. Journal of Chemical Information and Modeling, 2022, 62 (23), 5938-5951. DOI: 10.1021/acs.jcim.2c01073"""

## Prerequisites

Before building and running the project, ensure the following dependencies are installed:

- **CMake** (version 3.9 or higher)
- **GCC** or **Clang** (C++17 compatible compiler)
- **Git** (to clone the repository)
- **Eigen3** (see [getting started website](https://eigen.tuxfamily.org/dox/GettingStarted.html))

To install these dependencies on Ubuntu, you can use:

```bash
$ sudo apt update
$ sudo apt install cmake g++ python3 git
```

For other operating systems, refer to their respective package managers or installation guides.

## Installation
Clone the repository:
```bash
git clone https://github.com/your-repo/isim.git
cd isim
```

Create directory for build:

`$ mkdir build `

`$ cd build`

The project supports two build types: **Release** and **Debug**. By default, the build type is set to **Release**.

- To configure a **Release** build, use:

    `$ cmake -DCMAKE_BUILD_TYPE=Release .. `

    or simply:
    
    `$ cmake ..`

- To configure a **Debug** build, use:

    `$ cmake -DCMAKE_BUILD_TYPE=Debug ..`

Then compile the project with:

`$ make` 

After building, the binaries will be available in the `build` directory.

## Testing and Examples

Please refer to the `tests` directory for unit tests for the iSIM module. The tests ensure that the implementation is correct and that the results are consistent with expected outputs.

The `examples` directory contains example scripts demonstrating how to use the iSIM module for various applications, such as sampling, diversity selection, and clustering.

## Funding
Research contained in this package was supported by the National Institute of General Medical Sciences of the National Institutes  of  Health  under  award  number  R35GM150620.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
