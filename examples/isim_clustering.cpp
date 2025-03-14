#include <Eigen/Dense>
#include <iostream>
#include <fstream> 
#include "clustering.h"

// helper function to read fingerprints from a csv file
Eigen::ArrayXXf read_fps(std::string file_name){
    // Read the data from a csv file that contains fingerprints
    // a given row of the output contains the fingerprint of a molecule
    std::ifstream indata(file_name);
    if (!indata.is_open()) {
        throw std::runtime_error("Could not open file");
    }
    std::string line;
    std::vector<float> values; // Change to float
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stof(cell));
        }
        ++rows;
    }

    // Check if the file has only one column
    if (rows == values.size()) {
        Eigen::ArrayXXf array(rows, 1);
        for (uint i = 0; i < rows; ++i) {
            array(i, 0) = values[i];
        }
        return array;
    } else {
        Eigen::Map<const Eigen::ArrayXXf> array(values.data(), values.size() / rows, rows);
        return array.transpose();
    }
}

int main(){
    // Step 1: read the fingerprints from a csv file
    Eigen::ArrayXXf fps = read_fps("CHEMBL214_Ki_fps.csv").topRows(50); 
    // the file has to be in the same directory as the executable
    // the fingerprints must be binary, unless using the real functions.
    // the fingerprints must be in the format of n_fps x n_features
    // we are only using the first 50 rows for testing so that the example runs faster

    // Step 2.1: Initialize the HierarchicalClustering object
    // The HierarchicalClustering object is initialized with the fingerprints (fps).
    // This object will handle the clustering process.
    HierarchicalClustering hc = HierarchicalClustering(fps);

    // Step 2.2: Perform clustering
    // The clustering is performed using the "RR" similarity metric.
    // Other metrics like "JT" or "SM" can also be used depending on the requirement.
    hc.runClustering("RR");

    // Step 2.3: Get linkage matrix
    // The linkage matrix represents the hierarchical clustering structure.
    // Each row of the matrix contains information about the merged clusters.
    Eigen::ArrayXXf linkage_matrix = hc.getZ();

    // Step 2.4: Print the linkage matrix
    // The linkage matrix is printed to the console for verification or debugging purposes.
    std::cout << "Linkage matrix: " << std::endl;
    std::cout << linkage_matrix << std::endl;

    return 0;
}