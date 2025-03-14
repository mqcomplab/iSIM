#include "isim.h"
#include "hc_utils.h"
#include "comp.h"

class HierarchicalClustering{
    public:
        HierarchicalClustering(Eigen::ArrayXXf fps);
        void setTreeList(std::vector<HCTree> tree_lst);
        std::vector<HCTree> getTreeList();
        std::vector<int> maxIndices(std::string n_ary);
        HCTree runClustering(std::string n_ary);
        Eigen::ArrayXXf getZ();
        int getCurrentInd();

    private:
        std::vector<HCTree> tree_list; // Tree list that we want to perform HC on
        // Note: tree_list is modified during HC, in the end we get a list of only one tree (the final tree). 

        Eigen::ArrayXXf Z; // Z matrix 
        // column 1 of Z: idx of first cluster that is merged
        // column 2 of Z: idx of second cluster that is merged
        // column 3 of Z: number of the current merging step (feels a bit redundant)
        // column 4 of Z: number of molecules in the new/merged cluster

        int current_ind; // index available for new cluster
};

