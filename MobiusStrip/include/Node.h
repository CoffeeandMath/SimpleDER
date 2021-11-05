#ifndef NODE_H
#define NODE_H

#include <Eigen/Dense>
class Node {
public:
	Node();
	// Setting Node locations
	void setX(Eigen::Vector3d &Xset);
	Eigen::Vector3d getX();

	// Setting eta
	void seteta(double);
	double geteta();
	// Setting External force F
	void setFnode(Eigen::Vector3d F);
	Eigen::Vector3d getFnode();

	// Setting global node locations
	void setGlobalNodes(Eigen::VectorXi GN);
	Eigen::VectorXi getGlobalNodes();

	// Setting FreevsFixedNodes
	void FreeNodeLoc(bool);
	bool isXFree();
	void FreeEtaSet(bool);
	bool isEtaFree();
	void setMfl(Eigen::Matrix3d&);
	Eigen::Matrix3d getMfl();
	void setetaold(double);
	double getetaold();

protected:

private:
	Eigen::Vector3d X;
	Eigen::Vector3d Fnode;
	bool FreeNode;
	bool Freeeta;
	Eigen::VectorXi globalNodes;
	double eta;
	double etaold;
	Eigen::Matrix3d Mfl;

};

#endif // NODE_H
