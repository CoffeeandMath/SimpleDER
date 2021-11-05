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
#include <iostream>

class Material {
public:
	Material();
	void SetElementPointer(Element*);
	double getE();
	Eigen::Vector3d getdEdkappa();
	Eigen::Vector3d getdEdeta();
	Eigen::Matrix3d getddEdkappadeta();
	Eigen::Matrix3d getddEddeta();
	Eigen::Matrix3d getddEddkappa();
	void setkappa0(Eigen::Vector3d&);
	void Update(int);
	void setBMod(Eigen::VectorXd&);
	void setdS(double);
	void setw(double);
	void setKr(Eigen::Matrix3d&);
	void setnu(double);
	Eigen::Matrix3d getKr();
	double getetaprime();
	Eigen::Vector3d getdetaprime();
	Eigen::Matrix3d getddetaprime();
	double getw();

protected:

private:
	void calcE();
	void calcdEdkappa();
	void calcddEddkappa();
	void calcdEdeta();
	double etaprime();
	double etapprime();
	Eigen::Vector3d detapprime();
	double g(double);
	double dg(double);
	double ddg(double);
	Eigen::Vector3d detaprime();
	Eigen::Matrix3d ddetaprime();
	Eigen::Matrix3d cof(Eigen::Matrix3d&);

	Element *ElemP = NULL;

	double w;
	double nu;
	double window;
	double E;
	double dS;
	double dSinv;
	Eigen::Matrix3d Kr;
	Eigen::Matrix3d Qr;
	Eigen::Vector3d dEdkappa;
	Eigen::Matrix3d ddEddkappa;

	Eigen::Vector3d dEdeta;
	Eigen::Matrix3d ddEddeta;

	Eigen::Matrix3d ddEdkappadeta;

	Eigen::VectorXd BMod;
	Eigen::Vector3d kappa0;

};

#endif /* MATERIAL_H_ */
