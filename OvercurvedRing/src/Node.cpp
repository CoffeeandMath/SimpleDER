#include "Node.h"

Node::Node()
{
    X.setZero();

    Fnode.setZero();

    FreeNode = true;
    globalNodes.setZero();
}
void Node::setX(Eigen::Vector3d &Xset)
{
    X = Xset;
}
Eigen::Vector3d Node::getX()
{
    return X;
}

// Setting External force F
void Node::setFnode(Eigen::Vector3d F)
{
    Fnode = F;
}
Eigen::Vector3d Node::getFnode()
{
    return Fnode;
}

// Setting global node locations
void Node::setGlobalNodes(Eigen::VectorXi GN)
{
    globalNodes = GN;
}
Eigen::VectorXi Node::getGlobalNodes()
{
    return globalNodes;
}

// Setting FreevsFixedNodes
bool Node::isXFree(){
	return FreeNode;
}
void Node::FreeNodeLoc(bool fn)
{
    FreeNode = fn;
}
