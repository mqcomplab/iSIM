#include "isim.h"
#include "comp.h"

std::vector<int> getSortedIndices(const Eigen::ArrayXd& array);

std::vector<int> medoid_sampling(Eigen::ArrayXXf fps, std::string n_ary, double percentage);

std::vector<int> outlier_sampling(Eigen::ArrayXXf fps, std::string n_ary, double percentage);

std::vector<int> extremes_sampling(Eigen::ArrayXXf fps, std::string n_ary, double percentage);

std::vector<int> stratified_sampling(Eigen::ArrayXXf fps, std::string n_ary, double percentage, int strata);

std::vector<int> quota_sampling(Eigen::ArrayXXf fps, std::string n_ary, double percentage, int n_bins);