#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <sstream>
#include <gtest/gtest.h>

Eigen::ArrayXXf read_fps(std::string file_name);

Eigen::ArrayXXd read_compsim(std::string file_name);

class testIsim: public testing::Test{
    protected:
        testIsim(){
            //only use first 1000 rows for testing
            fps = read_fps("CHEMBL214_Ki_fps.csv").topRows(1000);
        }
        Eigen::ArrayXXf fps;
        double threshold = 1e-6; // threshold for comparing floating point numbers
};

void compare_arrays(Eigen::ArrayXf cpp_arr, Eigen::ArrayXf python_arr, double threshold);