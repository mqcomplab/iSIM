#include "div.h"
#include "test_utils.h"

TEST_F(testIsim, test_diversity){
    std::vector<int> idx_rr = diversity(fps, 20.0, "medoid", "RR", "isim");
    std::vector<int> idx_jt = diversity(fps, 60.0, "medoid", "JT", "isim");
    Eigen::ArrayXf idx_rr_python = read_fps("../python_results/diversity_RR_20percent_chembl214.csv");
    Eigen::ArrayXf idx_jt_python = read_fps("../python_results/diversity_JT_60percent_chembl214.csv");
    Eigen::ArrayXf idx_rr_eigen = Eigen::Map<Eigen::ArrayXi>(idx_rr.data(), idx_rr.size()).cast<float>();
    Eigen::ArrayXf idx_jt_eigen = Eigen::Map<Eigen::ArrayXi>(idx_jt.data(), idx_jt.size()).cast<float>();
    compare_arrays(idx_rr_eigen, idx_rr_python, threshold, true);
    compare_arrays(idx_jt_eigen, idx_jt_python, threshold, true);
}

TEST_F(testIsim, test_reverse_diversity){
    std::vector<int> idx_rr = reverse_diversity(fps, 20.0, "RR");
    std::vector<int> idx_jt = reverse_diversity(fps, 60.0, "JT");
    Eigen::ArrayXf idx_rr_python = read_fps("../python_results/reverse_diversity_RR_20percent_chembl214.csv");
    Eigen::ArrayXf idx_jt_python = read_fps("../python_results/reverse_diversity_JT_60percent_chembl214.csv");
    Eigen::ArrayXf idx_rr_eigen = Eigen::Map<Eigen::ArrayXi>(idx_rr.data(), idx_rr.size()).cast<float>();
    Eigen::ArrayXf idx_jt_eigen = Eigen::Map<Eigen::ArrayXi>(idx_jt.data(), idx_jt.size()).cast<float>();
    compare_arrays(idx_rr_eigen, idx_rr_python, threshold, true);
    compare_arrays(idx_jt_eigen, idx_jt_python, threshold, true);
}

