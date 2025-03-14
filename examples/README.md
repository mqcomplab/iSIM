# Examples

This folder contains example projects and scripts to demonstrate the usage and functionality of the iSIM framework. These examples are designed to help users get started quickly and understand how to integrate iSIM into their own projects.

## Contents

- **`isim_clustering.cpp`**: An example demonstrating the hierarchical clustering capabilities of the iSIM framework.
- **`isim_comp.cpp`**: A simple example with the `calculate_isim`, `pairwise_average`, and calculating the medoid or outlier for a set of points.
- **`isim_sampling.cpp`**: Contains examples for medoid, stratified, and quota sampling. 

For each example, update the CSV file name in the code to reference your own data file. These examples are self-contained and require minimal adjustments to function correctly.

## Compiling the examples

A simple `CMakeLists.txt` file is provided to compile the examples. To compile the example you are interested in, uncomment the corresponding `add_executable` line in the `CMakeLists.txt` file. All executables have the name `isim_example`. To build the examples, follow these steps:

1. Navigate to the `examples` directory:
   ```bash
   cd examples
   ```
2. Create a build directory:
   ```bash
    mkdir build
    cd build
    ```
3. Run CMake to configure the project:
    ```bash
    cmake ..
    ```
4. Build the examples:
    ```bash
    make
    ```
5. After building, you can run the examples:
    ```bash
    ./isim_example
    ```
### Example output for isim_sampling.cpp

```bash
$ ./isim_example
Medoid sampled indices: 
8 63 59 18 6 21 61 65 77 51 
Quota sampled indices: 
8 63 23 0 25 12 59 62 17 36 
Stratified sampled indices: 
8 63 60 28 84 85 19 44 96 45 
$
```

## Contributing

Feel free to contribute additional examples or improvements to the existing ones. Submit a pull request with your changes.

## License

This project is licensed under the [MIT License](../LICENSE).