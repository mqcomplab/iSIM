# Test Folder

This folder contains tests for the iSIM project. 
Currently the tests compare the C++ results on the CHEMBL214_KI.csv database to the python results. 
To download the dataset, go to: https://github.com/molML/MoleculeACE/tree/main/MoleculeACE/Data/benchmark_data/old 

## Dependencies

To execute the test cases, ensure that both **Eigen3** and **GoogleTest (GTest)** are installed as dependencies. For detailed installation guidance, visit the [Eigen3 Getting Started page](https://eigen.tuxfamily.org/dox/GettingStarted.html) and the [GoogleTest README](https://github.com/google/googletest/blob/main/googletest/README.md). Alternatively, you can use the *conda* package manager to install these dependencies within your environment.

## Running the tests

iSIM must be compiled in **Release mode** so that the flags for the test build and iSIM library match. 
In the test folder, run the following commands to build the test executable:

```
$ mkdir build
$ cd build
$ cmake ..
$ make 
```
In order to run the test executable, the CHEMBL214_KI_fps.csv file has to be in the same directory. After building the test executable, you can run the tests using the following command:
```
$ ./test
```

### Expected Output

When you run the tests, you should see output similar to the following:
``` 
[==========] Running 26 tests from 2 test suites.
[----------] Global test environment set-up.
[----------] 1 test from testClustering
[ RUN      ] testClustering.test_hc_clustering
Warning: No fingerprints provided.
Warning: No fingerprints provided.
Warning: No fingerprints provided.
[       OK ] testClustering.test_hc_clustering (410 ms)
[----------] 1 test from testClustering (410 ms total)

[----------] 25 tests from testIsim
[ RUN      ] testIsim.test_calculate_isim
[       OK ] testIsim.test_calculate_isim (370 ms)
[ RUN      ] testIsim.test_pairwise_average
[       OK ] testIsim.test_pairwise_average (820 ms)
[ RUN      ] testIsim.test_calculate_medoid
[       OK ] testIsim.test_calculate_medoid (393 ms)
```

\- The iSIM Development Team