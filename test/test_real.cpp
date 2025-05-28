#include "real.h"
#include "test_utils.h"

// these test cases are with binary FPS!!! 

TEST_F(testIsim, test_pairwise_average_real){
    double avg_rr = pairwise_average_real(fps, "RR");
    double avg_jt = pairwise_average_real(fps, "JT");
    double avg_sm = pairwise_average_real(fps, "SM");
    double diff_rr = std::abs(avg_rr - 0.20679641653372122);
    double diff_jt = std::abs(avg_jt - 0.33449004829932427);
    double diff_sm = std::abs(avg_sm - 0.5951895127549425);
    EXPECT_LE(diff_rr, threshold);
    EXPECT_LE(diff_jt, threshold);
    EXPECT_LE(diff_sm, threshold);
}

TEST_F(testIsim, test_calculate_isim_real){
    double isim_rr = calculate_isim_real(fps, "RR");
    double isim_jt = calculate_isim_real(fps, "JT");
    double isim_sm = calculate_isim_real(fps, "SM");
    double diff_rr = std::abs(isim_rr - 0.20679641653372122);
    double diff_jt = std::abs(isim_jt -0.3381198205187699);
    double diff_sm = std::abs(isim_sm - 0.5951895127549425);
    EXPECT_LE(diff_rr, threshold);
    EXPECT_LE(diff_jt, threshold);
    EXPECT_LE(diff_sm, threshold);
}

TEST_F(testIsim, test_calculate_compsim_real){
    Eigen::ArrayXd comp_sim_rr = calculate_comp_sim_real(fps, "RR");
    Eigen::ArrayXd comp_sim_jt = calculate_comp_sim_real(fps, "JT");
    Eigen::ArrayXd comp_sim_sm = calculate_comp_sim_real(fps, "SM");
    Eigen::ArrayXd comp_sim_rr_expected = read_compsim("../python_results/compsim_chembl214_1000_rows_RR_real.csv");
    Eigen::ArrayXd comp_sim_jt_expected = read_compsim("../python_results/compsim_chembl214_1000_rows_JT_real.csv");
    Eigen::ArrayXd comp_sim_sm_expected = read_compsim("../python_results/compsim_chembl214_1000_rows_SM_real.csv");
    EXPECT_TRUE(comp_sim_rr.isApprox(comp_sim_rr_expected, threshold));
    EXPECT_TRUE(comp_sim_jt.isApprox(comp_sim_jt_expected, threshold));
    EXPECT_TRUE(comp_sim_sm.isApprox(comp_sim_sm_expected, threshold));
}

TEST_F(testIsim, test_calculate_medoid_real){
    int medoid_rr = calculate_medoid_real(fps, "RR");
    int medoid_jt = calculate_medoid_real(fps, "JT");
    int medoid_sm = calculate_medoid_real(fps, "SM");
    EXPECT_EQ(medoid_rr, 8);
    EXPECT_EQ(medoid_jt, 8);
    EXPECT_EQ(medoid_sm, 81);
}

TEST_F(testIsim, test_calculate_outlier_real){
    int outlier_rr = calculate_outlier_real(fps, "RR");
    int outlier_jt = calculate_outlier_real(fps, "JT");
    int outlier_sm = calculate_outlier_real(fps, "SM");
    EXPECT_EQ(outlier_rr, 270);
    EXPECT_EQ(outlier_jt, 270);
    EXPECT_EQ(outlier_sm, 8);
}