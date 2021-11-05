/*
 * LocalConstructor.h
 *
 *  Created on: Jun 15, 2020
 *      Author: kevin
 */

#ifndef LOCALCONSTRUCTOR_H_
#define LOCALCONSTRUCTOR_H_
#include<Element.h>
#include<Material.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

class LocalConstructor {
public:
	LocalConstructor();

	double getE();
	Eigen::VectorXd getdE();
	Eigen::MatrixXd getddE();
	void Update(int);

    void setGlobalNodes(Eigen::VectorXi GN);
    Eigen::VectorXi getGlobalNodes();
    void setElementPointer(Element*);
    void setMaterialPointer(Material*);

protected:

private:

	void calcE();
	void calcdE();
	void calcddE();

	Element *ElemP = NULL;
	Material *MatP = NULL;

	double E;
	Eigen::VectorXd dE;
	Eigen::Matrix<double,11,11>  ddE;
	Eigen::VectorXi globalNodes;

};

#endif /* LOCALCONSTRUCTOR_H_ */
