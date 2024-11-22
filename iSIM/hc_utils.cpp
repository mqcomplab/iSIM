#include "hc_utils.h"

//Hierarchical clustering tree node methods

HCTreeNode::HCTreeNode(Eigen::ArrayXf fps, std::list<int> idx, int z_ind){
    setFps(fps);
    setIndices(idx);
    setNObjects(idx.size());
    setLeftPtr(NULL);
    setRightPtr(NULL);
    setIdx(z_ind);
};

Eigen::ArrayXf HCTreeNode::getFps(){
    return fingerprints;
}

void HCTreeNode::setFps(Eigen::ArrayXf fps){
    fingerprints = fps;
}

int HCTreeNode::getIdx(){
    return z_index;
}

void HCTreeNode::setIdx(int z_ind){
    z_index = z_ind;
}

HCTreeNode* HCTreeNode::getLeftPtr(){
    return leftPtr;
}

void HCTreeNode::setLeftPtr(HCTreeNode* input_left){
    leftPtr = input_left;
}

HCTreeNode* HCTreeNode::getRightPtr(){
    return rightPtr;
}

void HCTreeNode::setRightPtr(HCTreeNode* input_right){
    rightPtr = input_right;
}

int HCTreeNode::getNObjects(){
    return n_objects;
}

void HCTreeNode::setNObjects(int n_mol){
    n_objects = n_mol;
}

std::list<int> HCTreeNode::getIndices(){
    return obj_indices;
}

void HCTreeNode::setIndices(std::list<int> idx){
    obj_indices=idx;
}



//Hierarchical clustering tree class
HCTree::HCTree(){ 
    rootPtr = NULL;
}

HCTreeNode* HCTree::getRoot(){
    return rootPtr;
}

void HCTree::setRoot(HCTreeNode* input_root){
    rootPtr = input_root;
}

Eigen::ArrayXf HCTree::getRootFps(){
    return rootPtr->getFps();
}

int HCTree::getRootNObjects(){
    return rootPtr->getNObjects();
}

std::list<int> HCTree::getRootIndices(){
    return rootPtr->getIndices();
}

int HCTree::getRootIdx(){ //z_index
    return rootPtr->getIdx();
}

void HCTree::insertRoot(Eigen::ArrayXf fps, std::list<int> idx, int ind){
    rootPtr = new HCTreeNode(fps, idx, ind);
}

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