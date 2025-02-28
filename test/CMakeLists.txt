cmake_minimum_required(VERSION 3.7)

project(test)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "-g -O3 -march=native -fopenmp -Wall -Wextra -pipe -fPIC -flto=auto -fvisibility=hidden -pthread -mavx -mavx2 -msse4.2 -mtune=native")

find_package(Eigen3 REQUIRED)
find_package(GTest REQUIRED)

add_executable(test test_utils.cpp 
                test_clustering.cpp 
                test_comp.cpp 
                test_div.cpp
                test_real.cpp
                test_sampling.cpp
                test_main.cpp)
                
target_include_directories(test PRIVATE ../iSIM/ ${EIGEN3_INCLUDE_DIR})
target_link_directories(test PRIVATE ../build/)
target_link_libraries(test PRIVATE iSIM)
target_link_libraries(test PRIVATE iSIM GTest::GTest GTest::Main)