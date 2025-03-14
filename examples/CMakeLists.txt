cmake_minimum_required(VERSION 3.9)

project(Examples)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# uncomment to compile example for comp functions:
#add_executable(isim_example isim_comp.cpp)

# uncomment to compile example for clustering functions:
# add_executable(isim_example isim_clustering.cpp)

# uncomment to compile example for sampling functions:
add_executable(isim_example isim_sampling.cpp)

find_package(Eigen3 REQUIRED)

target_include_directories(isim_example
    PRIVATE 
        ${EIGEN3_INCLUDE_DIR}
        "../iSIM/" # relative path, can be changed to absolute if needed
)

target_link_directories(isim_example
    PRIVATE
        "../build/") # relative path, can be changed to absolute if needed

target_link_libraries(isim_example
    PRIVATE 
        iSIM
)

set(CMAKE_CXX_FLAGS "-O3 -march=native -Wall -Wextra -fopenmp")
