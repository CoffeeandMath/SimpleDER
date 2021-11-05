#include "Node.h"

Node::Node() {
	X.setZero();

	Fnode.setZero();
	Freeeta = true;
	FreeNode = true;
	globalNodes.setZero();
	Mfl.setZero();
	eta = 0.0;
	etaold = 0.0;
}
void Node::setX(Eigen::Vector3d &Xset) {
	X = Xset;
}
Eigen::Vector3d Node::getX() {
	return X;
}

// Setting External force F
void Node::setFnode(Eigen::Vector3d F) {
	Fnode = F;
}
Eigen::Vector3d Node::getFnode() {
	return Fnode;
}

// Setting global node locations
void Node::setGlobalNodes(Eigen::VectorXi GN) {
	globalNodes = GN;
}
Eigen::VectorXi Node::getGlobalNodes() {
	return globalNodes;
}

// Setting FreevsFixedNodes
bool Node::isXFree() {
	return FreeNode;
}
void Node::FreeNodeLoc(bool fn) {
	FreeNode = fn;
}
void Node::seteta(double et) {
	eta = et;
}
double Node::geteta() {
	return eta;
}
void Node::FreeEtaSet(bool fe) {
	Freeeta = fe;
}
bool Node::isEtaFree() {
	return Freeeta;
}
void Node::setMfl(Eigen::Matrix3d &Min) {
	Mfl = Min;

}
Eigen::Matrix3d Node::getMfl() {
	return Mfl;
}
void Node::setetaold(double d) {
	etaold = d;
}
double Node::getetaold() {
	return etaold;
}
