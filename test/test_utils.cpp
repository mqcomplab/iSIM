#include "test_utils.h"


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



Eigen::ArrayXXd read_compsim(std::string file_name) {
    // Read the data from a csv file that contains fingerprints
    // a given row of the output contains the fingerprint of a molecule
    std::ifstream indata(file_name);
    if (!indata.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::string line;
    std::vector<double> values; // Change to double
    uint rows = 0;

    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell)); // Change to stod
        }
        ++rows;
    }

    // Check if the file has only one column
    if (rows == values.size()) {
        Eigen::ArrayXXd array(rows, 1);
        for (uint i = 0; i < rows; ++i) {
            array(i, 0) = values[i];
        }
        return array;
    } else {
        Eigen::Map<const Eigen::ArrayXXd> array(values.data(), values.size() / rows, rows);
        return array.transpose();
    }
}

void compare_arrays(Eigen::ArrayXf cpp_arr, Eigen::ArrayXf python_arr, double threshold, bool sort){
    // sort arg is there so that we can compare if the chosen indices are the same regardless of the order
    // default is false
    ASSERT_EQ(cpp_arr.size(), python_arr.size());
    if (sort){
        std::sort(cpp_arr.begin(), cpp_arr.end());
        std::sort(python_arr.begin(), python_arr.end());
        bool equal = cpp_arr.isApprox(python_arr, threshold);
        EXPECT_TRUE(equal);
    }
    else{
        bool equal = cpp_arr.isApprox(python_arr, threshold);
        EXPECT_TRUE(equal);
    
        if (!equal){
            for (int i=0; i<cpp_arr.size(); i++){
                if (std::abs(cpp_arr(i)-python_arr(i)) > threshold){
                    std::cout << "index: " << i << ", c++ result: " << cpp_arr(i);
                    std::cout << ", python result: " << python_arr(i) << std::endl;
                }
            }
        }
    }
}