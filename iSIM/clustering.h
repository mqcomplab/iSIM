#include "isim.h"
#include "comp.h"

class HCTreeNode{
    // Class for the hierarchical clustering tree node
    // each node stores:
    //     fingerprints: the columwise sum of the fingerprints
    //     obj_indices: indices of the molecules in the cluster
    //     n_objects: number of molecules in the cluster
    // in addition to pointers to the left and right child nodes. 

    friend class Tree;

    private:
        int n_objects;
        Eigen::ArrayXf fingerprints;
        std::list<int> obj_indices;
        TreeNode* leftPtr;
        TreeNode* rightPtr;
        
    public:
        HCTreeNode(Eigen::ArrayXf fps, std::list<int> idx);

        Eigen::ArrayXf getFps();

        void setFps(Eigen::ArrayXf fps);

        HCTreeNode* getLeftPtr();

        void setLeftPtr(HCTreeNode* input_left);

        HCTreeNode* getRightPtr();

        void setRightPtr(HCTreeNode* input_right);

        void setNObjects(int n_mol);

        int getNObjects();

        void setIndices(std::list<int> idx);

        std::list<int> getIndices();
};

class HCTree{
    private:
        HCTreeNode* rootPtr;
    public:
        HCTree();

        HCTreeNode* getRoot();

        void setRoot(HCTreeNode* input_root);

        void insertRoot(Eigen::ArrayXf fps, std::list<int> idx);

        Eigen::ArrayXf getRootFps();

        int getRootNObjects();

        std::list<int> getRootIndices();


        void combineTrees(HCTree tree_left, HCTree tree_right);

        void printTree(HCTreeNode* root, int space);
};

