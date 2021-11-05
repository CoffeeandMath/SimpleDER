/*
 * Material.cpp
 *
 *  Created on: Jun 15, 2020
 *      Author: kevin
 */

#include "Material.h"

Material::Material() {
	// TODO Auto-generated constructor stub
	E = 0;
	//test
	dEdkappa.setZero();
	ddEddkappa.setZero();

	kappa0.setZero();
	BMod.setZero();
	li = 1.0;

}

double Material::Efunction(const Eigen::Vector3d & kappa){
	double Ev = 0.0;
	for (unsigned int i = 0; i < 3; i++)
	{
		Ev += 0.5*BMod(i)*pow(kappa(i) - kappa0(i),2.0);
	}
	return Ev;
}
Eigen::Vector3d Material::dEfunction(const Eigen::Vector3d & kappa){
	Eigen::Vector3d dEout;
	for (unsigned int i = 0; i < 3; i++)
		{
		dEout(i) = BMod(i)*(kappa(i) - kappa0(i));
		}
	return dEout;

}
Eigen::Matrix3d Material::ddEfunction(const Eigen::Vector3d & kappa){
	Eigen::Matrix3d ddEout;
	ddEout.setZero();

	for (unsigned int i = 0; i < 3; i++)
	{
		ddEout(i,i) = BMod(i);
	}

	return ddEout;
}

void Material::Update(int Order)
{
	calcE();
	if (Order > 0)
	{
		calcdEdkappa();
	}
	if (Order > 1)
	{
		calcddEddkappa();
	}
}

void Material::SetElementPointer(Element* Ein)
{
	ElemP = Ein;
}
void Material::setBMod(Eigen::Vector3d &Bm)
{
	BMod = Bm;
}
double Material::getE()
{
	return E;
}
Eigen::Vector3d Material::getdEdkappa()
{
	return dEdkappa;
}
Eigen::Matrix3d Material::getddEddkappa()
{
	return ddEddkappa;
}
void Material::setkappa0(Eigen::Vector3d &k0)
{
	kappa0 = k0;
}

void Material::calcE()
{
	E = li*Efunction((ElemP->getkappa())/li);
}
void Material::calcdEdkappa()
{
	dEdkappa = dEfunction((ElemP->getkappa())/li);

}
void Material::calcddEddkappa()
{
	ddEddkappa = ddEfunction((ElemP->getkappa())/li)/li;

}
void Material::setli(double & ltemp){
	li = ltemp;
}

