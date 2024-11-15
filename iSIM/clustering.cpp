#include "clustering.h"

HCTreeNode::HCTreeNode(Eigen::ArrayXf fps, std::list<int> idx){
    setFps(fps);
    setIndices(idx);
    setNObjects(idx.size());
    leftPtr = NULL;
    rightPtr = NULL;
};

Eigen::ArrayXf HCTreeNode::getFps(){
    return fingerprints;
}

void HCTreeNode::setFps(Eigen::ArrayXf fps){
    fingerprints = fps;
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

void HCTreeNode::setNObjects(int n_mol){
    n_objects = n_mol;
}

int HCTreeNode::getNObjects(){
    return n_objects;
}

void HCTreeNode::setIndices(std::list<int> idx){
    obj_indices=idx;
}

std::list<int> HCTreeNode::getIndices(){
    return obj_indices;
}



HCTree::HCTree(){
    root = NULL;
}

HCTreeNode* HCTree::getRoot(){
    return root;
}

void HCTree::setRoot(HCTreeNode* input_root){
    root = input_root;
}

void HCTree::insertRoot(Eigen::ArrayXf fps, std::list<int> idx){
    root = new TreeNode(fps, idx);
}

Eigen::ArrayXf HCTree::getRootFingerprints(){
    return root->getFps();
}

int HCTree::getRootNObjects(){
    return root->getNObjects();
}

std::list<int> HCTree::getRootIndices(){
    return root->getIndices();
}

void HCTree::combineTrees(HCTree tree_left, HCTree tree_right){
    Eigen::ArrayXf fps = tree_left.getRootFingerprints() + tree_right.getRootFingerprints();
    std::list<int> idx = tree_left.getRootIndices();
    idx.splice(idx.end(), tree_right.getRootIndices());
    idx.sort();
    TreeNode* new_root = new TreeNode(fps, idx);
    new_root->setLeft(tree_left.getRoot());
    new_root->setRight(tree_right.getRoot());
    
    setRoot(new_root);
}

void HCTree::printTree(HCTreeNode* root, int space){
    if (root==nullptr)
        return;
    space += COUNT;
    printTree(root->rightPtr, space);

    std::cout << std::endl;
    for (int i = COUNT; i < space; i++)
        std::cout << " ";
    std::list<int> idx = root->getIndices();
    std::cout << "(";
    for (int i : idx){
        std::cout << i << " ";
    }
    std::cout << ")\n";

    printTree(root->leftPtr, space);

}