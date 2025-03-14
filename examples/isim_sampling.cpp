#include <Eigen/Dense>
#include <iostream>
#include <fstream> 
#include "sampling.h"

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
    Eigen::ArrayXXf fps = read_fps("CHEMBL214_Ki_fps.csv").topRows(100); 
    // the file has to be in the same directory as the executable
    // the fingerprints must be binary, unless using the real functions.
    // the fingerprints must be in the format of n_fps x n_features
    // we are only using the first 1000 rows for testing so that the example runs faster

    // Step 2: Medoid sampling
    // sample 10% of the fingerprints using medoid sampling. Outlier and extreme sampling are similar. 
    std::vector<int> medoid_sampled = medoid_sampling(fps, "RR", 10);
    std::cout << "Medoid sampled indices: " << std::endl;
    for (int idx : medoid_sampled) {
        std::cout << idx << " ";
    }
    std::cout << std::endl;

    // Step 3: Quota sampling
    // sample 10% of the fingerprints using quota sampling with 5 bins
    std::vector<int> quota_sampled = quota_sampling(fps, "RR", 10, 5);
    std::cout << "Quota sampled indices: " << std::endl;
    for (int idx : quota_sampled) {
        std::cout << idx << " ";
    }
    std::cout << std::endl;

    // Step 5: Stratified sampling
    // sample 10% of the fingerprints using stratified sampling with 5 strata
    std::vector<int> stratified_sampled = stratified_sampling(fps, "RR", 10, 5);
    std::cout << "Stratified sampled indices: " << std::endl;
    for (int idx : stratified_sampled) {
        std::cout << idx << " ";
    }
    std::cout << std::endl;


    return 0;
}