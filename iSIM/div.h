#include <random>

#include "isim.h"
#include "comp.h"

std::vector<int> diversity(Eigen::ArrayXXf data, double percentage, std::string start, std::string n_ary, std::string method, int k=1);

std::vector<int> reverse_diversity(Eigen::ArrayXXf data, double percentage, std::string n_ary);