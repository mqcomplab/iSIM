#include "hc_utils.h"

/**
 * @brief Constructor for the HCTreeNode class.
 * @param fps A 1D Eigen Array representing the fingerprints for the node.
 * @param idx A list of integers representing the indices of the objects in this node.
 * @param z_ind An integer representing the z_index for this node. This is the index of the node in the linkage matrix Z.
 */
HCTreeNode::HCTreeNode(Eigen::ArrayXf fps, std::list<int> idx, int z_ind){
    setFps(fps);
    setIndices(idx);
    setNObjects(idx.size());
    setLeftPtr(NULL);
    setRightPtr(NULL);
    setIdx(z_ind);
};

/**
 * @brief Get the fingerprint for this node. Corresponds to the columnwise sum, if the node has multiple objects.
 * @return Eigen::ArrayXf The fingerprints for this node.
 */
Eigen::ArrayXf HCTreeNode::getFps(){
    return fingerprints;
}

/**
 * @brief Set the fingerprint for this node.
 * @param fps A 1D Eigen Array representing the fingerprints for the node.
 */
void HCTreeNode::setFps(Eigen::ArrayXf fps){
    fingerprints = fps;
}

/**
 * @brief Get the z_index for this node.
 * @return int The z_index for this node.
 */
int HCTreeNode::getIdx(){
    return z_index;
}

/**
 * @brief Set the z_index for this node.
 * @param z_ind An integer representing the z_index for this node.
 */
void HCTreeNode::setIdx(int z_ind){
    z_index = z_ind;
}

/**
 * @brief Get the left pointer for this node.
 * @return HCTreeNode* A pointer to the left child node.
 */
HCTreeNode* HCTreeNode::getLeftPtr(){
    return leftPtr;
}

/**
 * @brief Set the left pointer for this node.
 * @param input_left A pointer to the left child node.
 */
void HCTreeNode::setLeftPtr(HCTreeNode* input_left){
    leftPtr = input_left;
}

/**
 * @brief Get the right pointer for this node.
 * @return HCTreeNode* A pointer to the right child node.
 */
HCTreeNode* HCTreeNode::getRightPtr(){
    return rightPtr;
}

/**
 * @brief Set the right pointer for this node.
 * @param input_right A pointer to the right child node.
 */
void HCTreeNode::setRightPtr(HCTreeNode* input_right){
    rightPtr = input_right;
}

/**
 * @brief Get the number of objects in this node.
 * @return int The number of objects in this node.
 */
int HCTreeNode::getNObjects(){
    return n_objects;
}

/**
 * @brief Set the number of objects in this node.
 * @param n_mol An integer representing the number of objects in this node.
 */
void HCTreeNode::setNObjects(int n_mol){
    n_objects = n_mol;
}

/**
 * @brief Get the list of indices for the objects in this node.
 * @return std::list<int> A list of integers representing the indices of the objects in this node.
 */
std::list<int> HCTreeNode::getIndices(){
    return obj_indices;
}

/**
 * @brief Set the list of indices for the objects in this node.
 * @param idx A list of integers representing the indices of the objects in this node.
 */
void HCTreeNode::setIndices(std::list<int> idx){
    obj_indices=idx;
}



/**
 * @brief Constructor for the HCTree class.
 * 
 * @details
 * The HCTree class represents a hierarchical clustering tree. It contains a pointer to the root node
 * of the tree, which is an instance of the HCTreeNode class.
 * 
 */
HCTree::HCTree(){ 
    rootPtr = NULL;
}

/**
 * @brief Get the root node of the hierarchical clustering tree.
 * 
 * @return HCTreeNode* A pointer to the root node of the tree.
 */
HCTreeNode* HCTree::getRoot(){
    return rootPtr;
}

/**
 * @brief Set the root node of the hierarchical clustering tree.
 * 
 * @param input_root A pointer to the new root node of the tree.
 */
void HCTree::setRoot(HCTreeNode* input_root){
    rootPtr = input_root;
}

/**
 * @brief Get the fingerprints for the root node of the hierarchical clustering tree.
 * 
 * @return Eigen::ArrayXf The fingerprints for the root node of the tree.
 */
Eigen::ArrayXf HCTree::getRootFps(){
    return rootPtr->getFps();
}

/**
 * @brief Get the number of objects in the root node of the hierarchical clustering tree.
 * 
 * @return int The number of objects in the root node of the tree.
 */
int HCTree::getRootNObjects(){
    return rootPtr->getNObjects();
}

/**
 * @brief Get the list of indices for the objects in the root node of the hierarchical clustering tree.
 * 
 * @return std::list<int> A list of integers representing the indices of the objects in the root node of the tree.
 */
std::list<int> HCTree::getRootIndices(){
    return rootPtr->getIndices();
}

/**
 * @brief Get the z_index for the root node of the hierarchical clustering tree.
 * 
 * @return int The z_index for the root node of the tree.
 */
int HCTree::getRootIdx(){ //z_index
    return rootPtr->getIdx();
}

/**
 * @brief Insert a new root node into the hierarchical clustering tree.
 * 
 * @param fps A 1D Eigen Array representing the fingerprints for the new root node.
 * @param idx A list of integers representing the indices of the objects in the new root node.
 * @param ind An integer representing the z_index for the new root node.
 */
void HCTree::insertRoot(Eigen::ArrayXf fps, std::list<int> idx, int ind){
    rootPtr = new HCTreeNode(fps, idx, ind);
}

/**
 * @brief Combine two hierarchical clustering trees into one.
 * 
 * @details
 * This method combines two trees according to the hierarchical clustering algorithm. It calculates
 * the new fingerprints and indices for the combined tree, creates a new root node, and sets the
 * left and right pointers accordingly.
 * 
 * @param other_tree The other HCTree to be combined with this tree.
 * @param new_z_ind An integer representing the z_index for the new root node of the combined tree.
 */
void HCTree::combineTrees(HCTree other_tree, int new_z_ind){
    // This method combines two trees according, which is needed for hierarch. clust.
    
    // calculate new colsum for fingerprints, indices
    Eigen::ArrayXf fps = other_tree.getRootFps() + getRootFps();
    std::list<int> idx = other_tree.getRootIndices();
    idx.splice(idx.end(), getRootIndices());
    idx.sort();

    //create new root node
    HCTreeNode* new_root = new HCTreeNode(fps, idx, new_z_ind);

    // set left and right pointers so that z_left < z_right
    int z_other_tree = other_tree.getRootIdx();
    int z_this_tree = getRootIdx();
    if (z_this_tree < z_other_tree){
        new_root->setLeftPtr(getRoot());
        new_root->setRightPtr(other_tree.getRoot());
    }
    else{
        new_root->setLeftPtr(other_tree.getRoot());
        new_root->setRightPtr(getRoot());
    }
    // update root
    setRoot(new_root);
}