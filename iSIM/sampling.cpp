#include "sampling.h"
#include <vector>
#include <cmath>

//start only for testing:
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
// end only for testing

std::vector<int> getSortedIndices(const Eigen::ArrayXd& array) {
    std::vector<int> indices(array.size());
    for (int i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(), [&array](int a, int b) {
        return array[a] < array[b];
    });
    return indices;
}

std::vector<int> medoid_sampling(Eigen::ArrayXXf fps, std::string n_ary, double percentage){
    Eigen::ArrayXd comp_sims = calculate_comp_sim(fps, n_ary);
    std::vector<int> sorted_ind = getSortedIndices(comp_sims);
    int max_ind = sorted_ind.size()*percentage/100;
    std::vector<int> chosen_ind(max_ind);
    std::copy(sorted_ind.begin(), sorted_ind.begin() + max_ind, chosen_ind.begin());
    return chosen_ind;
}

std::vector<int> outlier_sampling(Eigen::ArrayXXf fps, std::string n_ary, double percentage){
    Eigen::ArrayXd comp_sims = calculate_comp_sim(fps, n_ary);
    std::vector<int> sorted_ind = getSortedIndices(comp_sims);
    int max_ind = sorted_ind.size()*percentage/100;
    std::vector<int> chosen_ind(max_ind);
    std::copy(sorted_ind.end()-max_ind, sorted_ind.end(), chosen_ind.begin());
    return chosen_ind;
}

std::vector<int> extremes_sampling(Eigen::ArrayXXf fps, std::string n_ary, double percentage){
    Eigen::ArrayXd comp_sims = calculate_comp_sim(fps, n_ary);
    std::vector<int> sorted_ind = getSortedIndices(comp_sims);
    int max_ind = sorted_ind.size()*percentage/200;
    std::vector<int> chosen_ind(max_ind);
    std::copy(sorted_ind.begin(), sorted_ind.begin()+ max_ind, chosen_ind.begin());
    std::vector<int> chosen_ind2(max_ind);
    std::copy(sorted_ind.end()-max_ind, sorted_ind.end(), chosen_ind2.begin());
    chosen_ind.insert(chosen_ind.end(), chosen_ind2.begin(), chosen_ind2.end());
    return chosen_ind;
}

std::vector<int> stratified_sampling(Eigen::ArrayXXf fps, std::string n_ary, double percentage, int strata){
    Eigen::ArrayXd comp_sims = calculate_comp_sim(fps, n_ary);
    std::vector<int> sorted_ind = getSortedIndices(comp_sims);
    int n_objects = sorted_ind.size();
    int n_stratum = std::floor(n_objects/strata);
    int rem_stratum = n_objects%strata;
    int n_sample = n_objects * percentage/100;
    if (n_sample < strata)
        throw std::invalid_argument("The number of objects to sample is smaller than the number of batches. Specify a higher percentage, or less batches.");
    int n_choose = std::floor(n_sample/strata);
    std::vector<int> chosen_ind;
    int rem = n_sample%strata;
    int stratum_idx = 0;
    for (int i = 0; i < strata; i++){
        // get elements from stratum
        if (i < rem){
            std::vector<int> chosen_idx_stratum(sorted_ind.begin()+stratum_idx, sorted_ind.begin() + stratum_idx + n_choose + 1);
            chosen_ind.insert(chosen_ind.end(), chosen_idx_stratum.begin(), chosen_idx_stratum.end());
        }
        else{
            std::vector<int> chosen_idx_stratum(sorted_ind.begin()+stratum_idx, sorted_ind.begin() + stratum_idx + n_choose);
            chosen_ind.insert(chosen_ind.end(), chosen_idx_stratum.begin(), chosen_idx_stratum.end());
        }
        // append chosen elements
        
        // update start index of stratum
        if (i < rem_stratum){
            stratum_idx += n_stratum+1;
        }
        else{
            stratum_idx += n_stratum;
        }
    }
    return chosen_ind;
}

std::vector<int> quota_sampling(Eigen::ArrayXXf fps, std::string n_ary, double percentage, int n_bins){
    Eigen::ArrayXd comp_sims = calculate_comp_sim(fps, n_ary);
    std::vector<int> sorted_ind = getSortedIndices(comp_sims);
    int n_objects = sorted_ind.size();
    int n_sample = n_objects * percentage/100;
    double min_comp_sim = comp_sims[sorted_ind[0]];
    double max_comp_sim = comp_sims[sorted_ind.back()];
    double comp_sim_diff = max_comp_sim - min_comp_sim;
    long double step = comp_sim_diff / n_bins;

    std::cout <<std::setprecision(std::numeric_limits<double>::digits10 + 1) << min_comp_sim << std::endl;
    std::cout << max_comp_sim << std::endl;
    std::cout << "step: " <<std::setprecision(std::numeric_limits<double>::digits10 + 1) << step << std::endl;
    std::vector<std::vector<int>> bins; // idx of molecules in a given bin
    int start_idx_bin = 0;
    for (int b_idx = 0; b_idx < n_bins+1; b_idx++){
        double bin_low = min_comp_sim + step*b_idx;
        double bin_high = min_comp_sim + step*(b_idx+1);
        std::vector<int> bin_b;
        for (int i = start_idx_bin; i < n_objects; i++){
            double comp_sim_i = comp_sims[sorted_ind[i]];
            if (comp_sim_i >= bin_low && comp_sim_i < bin_high){
                bin_b.push_back(sorted_ind[i]);
                // std::cout << comp_sims[sorted_ind[i]] << "     " << compsim_divide << std::endl;
            }
            else if(comp_sim_i >= bin_high){
                start_idx_bin = i;
                break;
            }
        }
        if (bin_b.size() != 0){
            bins.push_back(bin_b);
        }

}
    
    int n_bins_filled = bins.size();
    std::vector<int> sampled_idx;
    int i = 0;
    while (sampled_idx.size()<n_sample){
        for (int b=0; b < n_bins_filled; b++){
            if (bins[b].size()>i){
                sampled_idx.push_back(bins[b][i]);
                if(sampled_idx.size() >= n_sample){
                    break;
                }
            }
        }
        i++;
    }
    return sampled_idx;
}