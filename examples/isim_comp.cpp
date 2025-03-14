#include <Eigen/Dense>
#include <iostream>
#include <fstream> 
#include "comp.h"

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
    Eigen::ArrayXXf fps = read_fps("CHEMBL214_Ki_fps.csv").topRows(1000); 
    // the file has to be in the same directory as the executable
    // the fingerprints must be binary, unless using the real functions.
    // the fingerprints must be in the format of n_fps x n_features
    // we are only using the first 1000 rows for testing so that the example runs faster

    // Step 2: compute the instatnt similarity
    double isim = calculate_isim(fps, "RR"); 
    std::cout << "Instant similarity (RR): " << isim << std::endl;

    // Step 3: compute the pairwise average similarity
    double avg_isim = pairwise_average(fps, "RR");
    std::cout << "Pairwise average similarity (RR): " << avg_isim << std::endl;

    // Step 4: compute the medoid
    int medoid_index = calculate_medoid(fps, "RR");
    std::cout << "Medoid index (RR): " << medoid_index << std::endl;

    // Step 5: compute the outlier
    int outlier_index = calculate_outlier(fps, "RR");
    std::cout << "Outlier index (RR): " << outlier_index << std::endl;
    return 0;
}