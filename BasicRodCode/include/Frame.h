#ifndef FRAME_H
#define FRAME_H
#include <Node.h>
#include <quaternion.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>



class Frame
{
    public:
        Frame();
        void setNodesPointers(Node*,Node*);
        Node* getNode1Pointer();
        Node* getNode2Pointer();
        void setPhi(double);
        double getPhi();
        Eigen::Vector3d gett();
        Eigen::Vector3d getT();
        double geteps();
        void Update(int);
        void setD(const Quat<double>&);
        Quat<double> getD();
        Quat<double> getd();
        Quat<double> getptw();
        Quat<double> getppar();
        Eigen::Matrix3d getB();
        Eigen::Matrix<double, 3, 2> getdepsdx();
        Eigen::Tensor<double, 3> getDB();
        Eigen::Tensor<double, 4> getddepsdxdx();
        void setGlobalNodes(Eigen::VectorXi GN);
        Eigen::VectorXi getGlobalNodes();
        void FreePhi(bool fn);
        bool isPhiFree();

    protected:

    private:

   	void calct();
   	void calcT();
   	Quat<double> vtoQ(const Eigen::Vector3d &);
   	void calcptw();
   	void calcppar();
   	Eigen::Matrix3d skw(const Eigen::Vector3d &);

    Node *Node1Pointer = NULL;
    Node *Node2Pointer = NULL;

    double phi;
    double eps;
    Eigen::Vector3d t;
    Eigen::Vector3d T;
    Quat<double> D;
    Quat<double> ptw;
    Quat<double> ppar;
    Quat<double> d;
    Eigen::Matrix3d B;
    Eigen::Matrix<double,3,2> depsdx;
    Eigen::Tensor<double,3> DB;
    Eigen::Tensor<double,4> ddepsdxdx;
    Eigen::VectorXi globalNodes;
    bool FreePhibool;


};

#endif // FRAME_H
