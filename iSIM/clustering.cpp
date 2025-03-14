/**
 * @class HierarchicalClustering
 * @brief A class for performing hierarchical clustering on binary fingerprints.
 *
 * This class implements a hierarchical clustering algorithm using a bottom-up approach.
 * It takes a 2D Eigen Array of molecular fingerprints as input, where each row corresponds
 * to one molecule. The clustering process combines the most similar clusters iteratively
 * until a single cluster remains.
 *
 * @details
 * The class maintains a list of hierarchical clustering trees (HCTree) and a linkage matrix (Z).
 * The linkage matrix stores information about the clustering process, including the indices
 * of the merged clusters, the iteration at which they were merged, and the number of objects
 * in the resulting cluster.
 */

#include "clustering.h"

/**
 * @brief Constructor for the HierarchicalClustering class.
 * 
 * @param fps A 2D Eigen Array of molecular fingerprints, where each row corresponds to one molecule.
 */
HierarchicalClustering::HierarchicalClustering(Eigen::ArrayXXf fps){

    std::vector<HCTree> leaf_list= {};
    current_ind = 0;
    if (fps.rows() == 0){
        std::cerr << "Warning: No fingerprints provided." << std::endl;
    }
    for (int i = 0; i<fps.rows(); i++){
        HCTree tree = HCTree();
        std::list<int> idx = {i};
        tree.insertRoot(fps.row(i), idx, current_ind);
        leaf_list.push_back(tree);
        current_ind += 1;
    }
    setTreeList(leaf_list);
}

/**
 * @brief Get the linkage matrix Z.
 * 
 * @return Eigen::ArrayXXf The linkage matrix containing information about the clustering process.
 */
Eigen::ArrayXXf HierarchicalClustering::getZ(){
    return Z;
}

/**
 * @brief Set the tree list. 
 * 
 * @param tree_lst A vector of HCTree objects representing the current state of the clustering trees.
 */
void HierarchicalClustering::setTreeList(std::vector<HCTree> tree_lst){
    tree_list = tree_lst;
}

/**
 * @brief Get the current list of hierarchical clustering trees.
 * 
 * @return std::vector<HCTree> A vector of HCTree objects representing the current state of the clustering trees.
 */
std::vector<HCTree> HierarchicalClustering::getTreeList(){
    return tree_list;
}

/**
 * @brief Find the indices of the two most similar clusters based on their similarity score.
 * 
 * @param n_ary A string representing the n-ary type for similarity calculation. RR, JT, or SM. 
 * @return std::vector<int> A vector containing the indices of the two most similar clusters.
 */
std::vector<int> HierarchicalClustering::maxIndices(std::string n_ary){
    double max_sim = -3.08;
    int max1 = -1;
    int max2 = -1;
    std::vector<HCTree> current_tree_list = getTreeList();
    for (u_int i=0; i<current_tree_list.size(); i++){
        for(u_int j=i+1; j<current_tree_list.size(); j++){
            HCTree tree_i = current_tree_list[i];
            HCTree tree_j = current_tree_list[j];
            Eigen::ArrayXf fps_colsum = tree_i.getRootFps() + tree_j.getRootFps();
            int n = tree_i.getRootNObjects() + tree_j.getRootNObjects();
            double sim = calculate_isim(fps_colsum, n, n_ary); 
            if (sim>max_sim){
                max_sim = sim;
                max1 = i;
                max2 = j;
            }
        }
    }
    return {max1, max2};
}

/**
 * @brief Perform hierarchical clustering on the molecular fingerprints.
 * @note Currently there is no support for rerunning the clustering with a different n-ary type.
 * 
 * @param n_ary A string representing the n-ary type for similarity calculation. RR, JT, or SM. 
 * @return HCTree The final hierarchical clustering tree after all merges.
 */
HCTree HierarchicalClustering::runClustering(std::string n_ary){
    int N = getTreeList().size()-1;
    if (N<1){
        throw std::invalid_argument("Number of trees must be greater than 1.");
    }
    for (int i=0; i<N; i++){
        // combine most similar clusters
        std::vector<int> idx = maxIndices(n_ary);
        std::vector<HCTree> current_tree_list = getTreeList();
        HCTree tree_0 = current_tree_list[idx[0]];
        HCTree tree_1 = current_tree_list[idx[1]];
        int cluster_index1 = tree_0.getRootIdx();
        int cluster_index2 = tree_1.getRootIdx();
        tree_0.combineTrees(tree_1, current_ind);
        int n_comb = tree_0.getRootNObjects();

        // add data to Z matrix
        Eigen::ArrayXf zi(4);
        if (cluster_index1 < cluster_index2){
            zi << cluster_index1, cluster_index2, i+1, n_comb;
        }
        else{
            zi <<  cluster_index2, cluster_index1, i+1, n_comb;
        }
        Z.conservativeResize(Z.rows()+1, 4);
        Z.row(Z.rows()-1) = zi;
        current_ind += 1;

        // update tree list
        current_tree_list[idx[0]] = tree_0;
        current_tree_list.erase(current_tree_list.begin() + idx[1]);
        setTreeList(current_tree_list);
    }
    HCTree final_tree = getTreeList()[0];
    return final_tree;
}