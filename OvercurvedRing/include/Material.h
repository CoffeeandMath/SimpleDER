/*
 * Material.h
 *
 *  Created on: Jun 15, 2020
 *      Author: kevin
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Element.h>
#include <math.h>

class Material {
public:
	Material();
	void SetElementPointer(Element*);
	double getE();
	Eigen::Vector3d getdEdkappa();
	Eigen::Matrix3d getddEddkappa();
	void setkappa0(Eigen::Vector3d &);
	void Update(int);
	void setBMod(Eigen::Vector3d &);


protected:

private:
	void calcE();
	void calcdEdkappa();
	void calcddEddkappa();


    Element *ElemP = NULL;

	double E;
	Eigen::Vector3d dEdkappa;
	Eigen::Matrix3d ddEddkappa;

	Eigen::Vector3d BMod;
	Eigen::Vector3d kappa0;

};

#endif /* MATERIAL_H_ */
