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
	dEdkappa.setZero();
	ddEddkappa.setZero();

	kappa0.setZero();
	BMod.setZero();

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
	E = 0.0;
	for (unsigned int i = 0; i < 3; i++)
	{
		E += 0.5*BMod(i)*pow((ElemP->getkappa())(i) - kappa0(i),2.0);
	}
}
void Material::calcdEdkappa()
{
	for (unsigned int i = 0; i < 3; i++)
	{
		dEdkappa(i) = BMod(i)*((ElemP->getkappa())(i) - kappa0(i));
	}

}
void Material::calcddEddkappa()
{
	//ddEddkappa.setZero();
	for (unsigned int i = 0; i < 3; i++)
	{
		ddEddkappa(i,i) = BMod(i);
	}

}

