#include "sampling.h"
#include "test_utils.h"

TEST_F(testIsim, test_medoid_sampling_RR){
    std::vector<int> idx_rr = medoid_sampling(fps, "RR", 20.0);
    Eigen::ArrayXf idx_python_rr = read_fps("../python_results/medoid_sampling_RR_20percent_chembl214.csv");
    Eigen::ArrayXf idx_eigen_rr = Eigen::Map<Eigen::ArrayXi>(idx_rr.data(), idx_rr.size()).cast<float>();
    compare_arrays(idx_python_rr, idx_eigen_rr, threshold, true);
}

TEST_F(testIsim, test_medoid_sampling_JT){
    std::vector<int> idx_jt = medoid_sampling(fps, "JT", 35.0);
    Eigen::ArrayXf idx_python_jt = read_fps("../python_results/medoid_sampling_JT_35percent_chembl214.csv");
    Eigen::ArrayXf idx_eigen_jt = Eigen::Map<Eigen::ArrayXi>(idx_jt.data(), idx_jt.size()).cast<float>();
    compare_arrays(idx_python_jt, idx_eigen_jt, threshold, true);
}

TEST_F(testIsim, test_medoid_sampling_SM){
    std::vector<int> idx_sm = medoid_sampling(fps, "SM", 50.0);
    Eigen::ArrayXf idx_python_sm = read_fps("../python_results/medoid_sampling_SM_50percent_chembl214.csv");
    Eigen::ArrayXf idx_eigen_sm = Eigen::Map<Eigen::ArrayXi>(idx_sm.data(), idx_sm.size()).cast<float>();
    compare_arrays(idx_python_sm, idx_eigen_sm, threshold,true);
}

TEST_F(testIsim, test_outlier_sampling_RR){
    std::vector<int> idx_rr = outlier_sampling(fps, "RR", 35.0);
    Eigen::ArrayXf idx_python_rr = read_fps("../python_results/outlier_sampling_RR_35percent_chembl214.csv");
    Eigen::ArrayXf idx_eigen_rr = Eigen::Map<Eigen::ArrayXi>(idx_rr.data(), idx_rr.size()).cast<float>();
    compare_arrays(idx_python_rr, idx_eigen_rr, threshold,true);
}

TEST_F(testIsim, test_outlier_sampling_JT){
    std::vector<int> idx_jt = outlier_sampling(fps, "JT", 50.0);
    Eigen::ArrayXf idx_python_jt = read_fps("../python_results/outlier_sampling_JT_50percent_chembl214.csv");
    Eigen::ArrayXf idx_eigen_jt = Eigen::Map<Eigen::ArrayXi>(idx_jt.data(), idx_jt.size()).cast<float>();
    compare_arrays(idx_python_jt, idx_eigen_jt, threshold,true);
}

TEST_F(testIsim, test_outlier_sampling_SM){
    std::vector<int> idx_sm = outlier_sampling(fps, "SM", 20.0);
    Eigen::ArrayXf idx_python_sm = read_fps("../python_results/outlier_sampling_SM_20percent_chembl214.csv");
    Eigen::ArrayXf idx_eigen_sm = Eigen::Map<Eigen::ArrayXi>(idx_sm.data(), idx_sm.size()).cast<float>();
    compare_arrays(idx_python_sm, idx_eigen_sm, threshold,true);
}

TEST_F(testIsim, test_extremes_sampling_RR){
    std::vector<int> idx_rr = extremes_sampling(fps, "RR", 50.0);
    Eigen::ArrayXf idx_python_rr = read_fps("../python_results/extremes_sampling_RR_50percent_chembl214.csv");
    Eigen::ArrayXf idx_eigen_rr = Eigen::Map<Eigen::ArrayXi>(idx_rr.data(), idx_rr.size()).cast<float>();
    compare_arrays(idx_python_rr, idx_eigen_rr, threshold,true);
}

TEST_F(testIsim, test_extremes_sampling_JT){
    std::vector<int> idx_jt = extremes_sampling(fps, "JT", 20.0);
    Eigen::ArrayXf idx_python_jt = read_fps("../python_results/extremes_sampling_JT_20percent_chembl214.csv");
    Eigen::ArrayXf idx_eigen_jt = Eigen::Map<Eigen::ArrayXi>(idx_jt.data(), idx_jt.size()).cast<float>();
    compare_arrays(idx_python_jt, idx_eigen_jt, threshold, true);
}

TEST_F(testIsim, test_extremes_sampling_SM){
    std::vector<int> idx_sm = extremes_sampling(fps, "SM", 35.0);
    Eigen::ArrayXf idx_python_sm = read_fps("../python_results/extremes_sampling_SM_35percent_chembl214.csv");
    Eigen::ArrayXf idx_eigen_sm = Eigen::Map<Eigen::ArrayXi>(idx_sm.data(), idx_sm.size()).cast<float>();
    compare_arrays(idx_python_sm, idx_eigen_sm, threshold, true);
}

// only using RR for stratified and quota sampling

TEST_F(testIsim, test_stratified_sampling_RR){
    std::vector<int> idx_strat = stratified_sampling(fps, "RR", 40.0, 15);
    Eigen::ArrayXf idx_python_strat = read_fps("../python_results/stratified_sampling_RR_40percent_15strata_chembl214.csv");
    Eigen::ArrayXf idx_eigen_strat = Eigen::Map<Eigen::ArrayXi>(idx_strat.data(), idx_strat.size()).cast<float>();
    compare_arrays(idx_python_strat, idx_eigen_strat, threshold, true);
}

TEST_F(testIsim, test_quota_sampling_RR){
    std::vector<int> idx_quota = quota_sampling(fps, "RR", 40.0, 15);
    Eigen::ArrayXf idx_python_quota = read_fps("../python_results/quota_sampling_RR_40percent_15bins_chembl214.csv");
    Eigen::ArrayXf idx_eigen_quota = Eigen::Map<Eigen::ArrayXi>(idx_quota.data(), idx_quota.size()).cast<float>();
    compare_arrays(idx_python_quota, idx_eigen_quota, threshold, true);
}