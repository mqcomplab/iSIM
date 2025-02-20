#include "clustering.h"
#include "test_utils.h"

class testClustering: public testing::Test{
    protected:
        testClustering(){
            // only use first 50 rows for testing
            fps = read_fps("CHEMBL214_Ki_fps.csv").topRows(50);
            hc_rr = HierarchicalClustering(fps);
            hc_jt = HierarchicalClustering(fps);
            hc_sm = HierarchicalClustering(fps);
        }
        Eigen::ArrayXXf fps;
        double threshold = 1e-6; // threshold for comparing floating point numbers
        HierarchicalClustering hc_rr  = HierarchicalClustering(fps);
        HierarchicalClustering hc_jt  = HierarchicalClustering(fps);
        HierarchicalClustering hc_sm  = HierarchicalClustering(fps);
};

TEST_F(testClustering, test_hc_clustering){
    hc_rr.runClustering("RR");
    Eigen::ArrayXXf Z_rr = hc_rr.getZ();
    hc_jt.runClustering("JT");
    Eigen::ArrayXXf Z_jt = hc_jt.getZ();
    hc_sm.runClustering("SM");
    Eigen::ArrayXXf Z_sm = hc_sm.getZ();
    Eigen::ArrayXXf Z_rr_python = read_fps("../python_results/hc_clustering_50_fps_RR_chembl214.csv");
    Eigen::ArrayXXf Z_jt_python = read_fps("../python_results/hc_clustering_50_fps_JT_chembl214.csv");
    Eigen::ArrayXXf Z_sm_python = read_fps("../python_results/hc_clustering_50_fps_SM_chembl214.csv");
    ASSERT_EQ(Z_rr.rows(), Z_rr_python.rows());
    EXPECT_TRUE(Z_rr.isApprox(Z_rr_python, threshold));
    ASSERT_EQ(Z_jt.rows(), Z_jt_python.rows());
    EXPECT_TRUE(Z_jt.isApprox(Z_jt_python, threshold));
    ASSERT_EQ(Z_sm.rows(), Z_sm_python.rows());
    EXPECT_TRUE(Z_sm.isApprox(Z_sm_python, threshold));
}