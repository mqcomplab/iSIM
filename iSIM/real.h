#include "isim.h"

double pairwise_average_real(const Eigen::ArrayXXf fingerprints, std::string n_ary);

double calculate_isim_real(const Eigen::ArrayXXf fingerprints, std::string n_ary);

Eigen::ArrayXf calculate_comp_sim_real(const Eigen::ArrayXXf fingerprints, std::string n_ary);

int calculate_medoid_real(const Eigen::ArrayXXf data, const std::string n_ary);

int calculate_outlier_real(const Eigen::ArrayXXf data, const std::string n_ary);