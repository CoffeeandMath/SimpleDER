#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <Node.h>
#include <quaternion.h>
#include <Frame.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>



class Element {
public:
	Element();
	void setFramesPointers(Frame*,Frame*);
	void setNodesPointers(Node*,Node*,Node*);
	void Update(int);
	Quat<double> getq();
	Eigen::Vector3d getkappa();
    Eigen::Matrix<double,3,2> getdkappadphi();
    Eigen::Tensor<double,3> getdkappadx();
    Eigen::Tensor<double,3> getddkappadphidphi();
    Eigen::Tensor<double,4> getddkappadxdphi();
    Eigen::Tensor<double,5> getddkappadxdx();

    void setGlobalNodes(Eigen::VectorXi GN);
    Eigen::VectorXi getGlobalNodes();

protected:

private:

    void calcdkappadx();
    void calcdkappadphi();
    void calcddkappadphidphi();
    void calcddkappadxdphi();
    void calcddkappadxdx();

    Frame *Frame1P = NULL;
    Frame *Frame2P = NULL;

    Node *Node1P = NULL;
    Node *Node2P = NULL;
    Node *Node3P = NULL;

    Eigen::VectorXi globalNodes;

    Quat<double> q;
    Eigen::Vector3d kappa;
    Eigen::Matrix<double,3,2> dkappadphi;
    Eigen::Tensor<double,3> dkappadx;
    Eigen::Tensor<double,3> ddkappadphidphi;
    Eigen::Tensor<double,4> ddkappadxdphi;
    Eigen::Tensor<double,5> ddkappadxdx;

    Quat<double> qe3;



};

#endif /* ELEMENT_H_ */
