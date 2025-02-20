#include "comp.h"
#include "test_utils.h"



TEST_F(testIsim, test_calculate_isim){
    // python results are for databae CHEMBL214_Ki_fps.csv with 1000 rows
    double isim_rr = calculate_isim(fps, "RR");
    double diff_rr = std::abs(isim_rr - 0.20679641653372122); 
    EXPECT_LE(diff_rr, threshold); 
    double isim_jt = calculate_isim(fps, "JT");
    double diff_jt = std::abs(isim_jt - 0.3381198205187699);
    EXPECT_LE(diff_jt, threshold);
    double isim_sm = calculate_isim(fps, "SM");
    double diff_sm = std::abs(isim_sm - 0.5951895127549425);
    EXPECT_LE(diff_sm, threshold);
}

TEST_F(testIsim, test_pairwise_average){
    double avg_rr = pairwise_average(fps, "RR");
    double diff_rr = std::abs(avg_rr - 0.20679641653372122);
    EXPECT_LE(diff_rr, threshold);
    double avg_jt = pairwise_average(fps, "JT");
    double diff_jt = std::abs(avg_jt - 0.33449004829932427); 
    EXPECT_LE(diff_jt, threshold);
    double avg_sm = pairwise_average(fps, "SM");
    double diff_sm = std::abs(avg_sm - 0.5951895127549425);
    EXPECT_LE(diff_sm, threshold);
}

TEST_F(testIsim, test_calculate_medoid){
    int medoid_rr = calculate_medoid(fps, "RR");
    EXPECT_EQ(medoid_rr, 8);
    int medoid_jt = calculate_medoid(fps, "JT");
    EXPECT_EQ(medoid_jt, 8);
    int medoid_sm = calculate_medoid(fps, "SM");
    EXPECT_EQ(medoid_sm, 81);
}

TEST_F(testIsim, test_calculate_outlier){
    int outlier_rr = calculate_outlier(fps, "RR");
    EXPECT_EQ(outlier_rr, 270);
    int outlier_jt = calculate_outlier(fps, "JT");
    EXPECT_EQ(outlier_jt, 270);
    int outlier_sm = calculate_outlier(fps, "SM");
    EXPECT_EQ(outlier_sm, 8);
}

TEST_F(testIsim, test_calculate_counters){
    // testing relative differences, since the numbers here are very large
    int k = 3;
    std::map<std::string, double> counters = calculate_counters(fps, k);
    double rel_diff_a = std::abs(counters["a"] - 87219.84572536658)/87219.84572536658;
    EXPECT_LE(rel_diff_a, threshold);
    double rel_diff_d = std::abs(counters["d"] - 112447.45934277849)/112447.45934277849;
    EXPECT_LE(rel_diff_d, threshold);
    double rel_diff_total_dis = std::abs(counters["total_dis"] -119054.09706176948)/119054.09706176948;
    EXPECT_LE(rel_diff_total_dis, threshold);
    double rel_diff_total_sim = std::abs(counters["total_sim"] - 199667.30506814507)/199667.30506814507;
    EXPECT_LE(rel_diff_total_sim, threshold);
    double rel_diff_p = std::abs(counters["p"] - 318721.40212991455)/318721.40212991455;
    EXPECT_LE(rel_diff_p, threshold);
}

TEST_F(testIsim, test_calculate_compsim){
    Eigen::ArrayXd comp_sim_rr = calculate_comp_sim(fps, "RR");
    Eigen::ArrayXd comp_sim_rr_python = read_compsim("../python_results/compsim_chembl214_1000_rows_RR.csv");
    Eigen::ArrayXd comp_sim_jt = calculate_comp_sim(fps, "JT");
    Eigen::ArrayXd comp_sim_jt_python = read_compsim("../python_results/compsim_chembl214_1000_rows_JT.csv");
    Eigen::ArrayXd comp_sim_sm = calculate_comp_sim(fps, "SM");
    Eigen::ArrayXd comp_sim_sm_python = read_compsim("../python_results/compsim_chembl214_1000_rows_SM.csv");

    EXPECT_TRUE(comp_sim_rr.isApprox(comp_sim_rr_python, threshold));
    EXPECT_TRUE(comp_sim_jt.isApprox(comp_sim_jt_python, threshold));
    EXPECT_TRUE(comp_sim_sm.isApprox(comp_sim_sm_python, threshold));
}

TEST_F(testIsim, test_gen_sim_dict){
    std::map<std::string, double> sim_dict = gen_sim_dict(fps, 3);
    double diff_ac = std::abs(sim_dict["AC"] - 0.581393207617719);
    EXPECT_LE(diff_ac, threshold);
    double diff_bub = std::abs(sim_dict["BUB"] - 0.6100518722868215);
    EXPECT_LE(diff_bub, threshold);
    double diff_fai = std::abs(sim_dict["Fai"] - 0.4500594388646877);
    EXPECT_LE(diff_fai, threshold);
    double diff_gle = std::abs(sim_dict["Gle"] - 0.5943556500286961);
    EXPECT_LE(diff_gle, threshold);
    double diff_ja = std::abs(sim_dict["Ja"] - 0.6872870148186374);
    EXPECT_LE(diff_ja, threshold);
    double diff_jt = std::abs(sim_dict["JT"] - 0.4228350151592967);
    EXPECT_LE(diff_jt, threshold);
    double diff_rt = std::abs(sim_dict["RT"] - 0.4560952027621785);
    EXPECT_LE(diff_rt, threshold);
    double diff_rr = std::abs(sim_dict["RR"] - 0.2736554405901326);
    EXPECT_LE(diff_rr, threshold);
    double diff_sm = std::abs(sim_dict["SM"] - 0.6264634371392429);
    EXPECT_LE(diff_sm, threshold);
    double diff_ss1 = std::abs(sim_dict["SS1"] - 0.2680981503035359);
    EXPECT_LE(diff_ss1, threshold);
    double diff_ss2 = std::abs(sim_dict["SS2"] - 0.7703381740214438);
    EXPECT_LE(diff_ss2, threshold);
}