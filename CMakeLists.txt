cmake_minimum_required(VERSION 3.7)

project(iSIM 
    LANGUAGES CXX 
    VERSION 0.1.0
    DESCRIPTION "Module to perform instant similarity (iSIM), which allows for comparison of multiple molecules at the same time, yielding the same value as the average pairwise comparisons of other methods."
    HOMEPAGE_URL "https://github.com/mqcomplab/iSIM"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3 -march=native -fopenmp -Wall -Wextra -pipe -fPIC -flto=auto -fvisibility=hidden -pthread -mavx -mavx2 -msse4.2 -mtune=native")

add_library(iSIM  STATIC
    iSIM/clustering.cpp
    iSIM/comp.cpp
    iSIM/div.cpp
    iSIM/hc_utils.cpp
    iSIM/real.cpp 
    iSIM/sampling.cpp 
)

find_package(Eigen3 REQUIRED)
target_include_directories(iSIM PRIVATE ${EIGEN3_INCLUDE_DIR})
