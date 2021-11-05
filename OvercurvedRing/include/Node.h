#ifndef NODE_H
#define NODE_H

#include <Eigen/Dense>
class Node
{
public:
    Node();
    // Setting Node locations
    void setX(Eigen::Vector3d &Xset);
    Eigen::Vector3d getX();

    // Setting External force F
    void setFnode(Eigen::Vector3d F);
    Eigen::Vector3d getFnode();

    // Setting global node locations
    void setGlobalNodes(Eigen::VectorXi GN);
    Eigen::VectorXi getGlobalNodes();

    // Setting FreevsFixedNodes
    void FreeNodeLoc(bool fn);
    bool isXFree();



protected:

private:
    Eigen::Vector3d X;
    Eigen::Vector3d Fnode;
    bool FreeNode;
    Eigen::VectorXi globalNodes;


};

#endif // NODE_H
