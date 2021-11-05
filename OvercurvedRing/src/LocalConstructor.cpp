/*
 * LocalConstructor.cpp
 *
 *  Created on: Jun 15, 2020
 *      Author: kevin
 */

#include "LocalConstructor.h"

LocalConstructor::LocalConstructor() {
	// TODO Auto-generated constructor stub
	E = 0.0;
	dE.resize(11);
	dE.setZero();
	ddE.setZero();
	globalNodes.resize(11);
	globalNodes.setZero();
}

double LocalConstructor::getE()
{
	return E;
}

Eigen::VectorXd LocalConstructor::getdE()
{
	return dE;
}

Eigen::MatrixXd LocalConstructor::getddE()
{
	return ddE;
}

void LocalConstructor::Update(int Order)
{
	calcE();
	if (Order > 0)
	{
		calcdE();
	}
	if (Order > 1)
	{
		calcddE();
	}
}

void LocalConstructor::calcE()
{
	E = (MatP->getE());
}
void LocalConstructor::calcdE()
{
	Eigen::VectorXd dEdphi = (ElemP->getdkappadphi()).transpose()*(MatP->getdEdkappa());
	Eigen::Vector3d dEdxim; dEdxim.setZero();
	Eigen::Vector3d dEdxi; dEdxi.setZero();
	Eigen::Vector3d dEdxip; dEdxip.setZero();

	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			dEdxim(i) += (ElemP->getdkappadx())(j,i,0)*(MatP->getdEdkappa())(j);
			dEdxi(i) += (ElemP->getdkappadx())(j,i,1)*(MatP->getdEdkappa())(j);
			dEdxip(i) += (ElemP->getdkappadx())(j,i,2)*(MatP->getdEdkappa())(j);
		}

	}
	//dE << dEdxim(0), dEdxim(1), dEdxim(2), dEdphi(0), dEdxi(0), dEdxi(1), dEdxi(2), dEdphi(1), dEdxip(0), dEdxip(1), dEdxip(2);
	dE << dEdxim, dEdphi(0), dEdxi, dEdphi(1), dEdxip;
}
void LocalConstructor::calcddE()
{
	ddE.setZero();

	//dxdx
	//Second derivative of E(kappa) term
	int n1off = 0;
	int n2off = 4;
	int n3off = 8;
	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			for (unsigned int k = 0; k < 3; k++)
			{
				for (unsigned int l = 0; l < 3; l++)
				{
					ddE(i+n1off,l+n1off) += (ElemP->getdkappadx())(j,i,0)*(MatP->getddEddkappa())(j,k)*(ElemP->getdkappadx())(k,l,0);
					ddE(i+n1off,l+n2off) += (ElemP->getdkappadx())(j,i,0)*(MatP->getddEddkappa())(j,k)*(ElemP->getdkappadx())(k,l,1);
					ddE(i+n1off,l+n3off) += (ElemP->getdkappadx())(j,i,0)*(MatP->getddEddkappa())(j,k)*(ElemP->getdkappadx())(k,l,2);
					ddE(i+n2off,l+n2off) += (ElemP->getdkappadx())(j,i,1)*(MatP->getddEddkappa())(j,k)*(ElemP->getdkappadx())(k,l,1);
					ddE(i+n2off,l+n3off) += (ElemP->getdkappadx())(j,i,1)*(MatP->getddEddkappa())(j,k)*(ElemP->getdkappadx())(k,l,2);
					ddE(i+n3off,l+n3off) += (ElemP->getdkappadx())(j,i,2)*(MatP->getddEddkappa())(j,k)*(ElemP->getdkappadx())(k,l,2);

				}
			}
		}
	}

	//Second derivative of kappa term
	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			for (unsigned int l = 0; l < 3; l++)
			{

					ddE(i+n1off,l+n1off) += (ElemP->getddkappadxdx())(j,i,l,0,0)*(MatP->getdEdkappa())(j);
					ddE(i+n1off,l+n2off) += (ElemP->getddkappadxdx())(j,i,l,0,1)*(MatP->getdEdkappa())(j);
					ddE(i+n1off,l+n3off) += (ElemP->getddkappadxdx())(j,i,l,0,2)*(MatP->getdEdkappa())(j);
					ddE(i+n2off,l+n2off) += (ElemP->getddkappadxdx())(j,i,l,1,1)*(MatP->getdEdkappa())(j);
					ddE(i+n2off,l+n3off) += (ElemP->getddkappadxdx())(j,i,l,1,2)*(MatP->getdEdkappa())(j);
					ddE(i+n3off,l+n3off) += (ElemP->getddkappadxdx())(j,i,l,2,2)*(MatP->getdEdkappa())(j);

			}
		}
	}
	//Symmetries
	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			ddE(i+n2off,j+n1off) = ddE(j+n1off,i+n2off);
			ddE(i+n3off,j+n1off) = ddE(j+n1off,i+n3off);
			ddE(i+n3off,j+n2off) = ddE(j+n2off,i+n3off);
		}
	}

	//dxdphi
	//Second derivative of E(kappa) term
	Eigen::Vector2i philocs;
	philocs << 3,7;

	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			for (unsigned int k = 0; k < 3; k++)
			{
				for (unsigned int l = 0; l < 2; l++)
				{
					ddE(i+n1off,philocs(l)) += (ElemP->getdkappadx())(j,i,0) * (MatP->getddEddkappa())(j,k) * (ElemP->getdkappadphi())(k,l);
					ddE(i+n2off,philocs(l)) += (ElemP->getdkappadx())(j,i,1) * (MatP->getddEddkappa())(j,k) * (ElemP->getdkappadphi())(k,l);
					ddE(i+n3off,philocs(l)) += (ElemP->getdkappadx())(j,i,2) * (MatP->getddEddkappa())(j,k) * (ElemP->getdkappadphi())(k,l);
				}
			}

		}
	}

	//Second derivative of kappa term
	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			for (unsigned int l = 0; l < 2; l++)
			{
				ddE(i+n1off,philocs(l)) += (MatP->getdEdkappa())(j) * (ElemP->getddkappadxdphi())(j,i,0,l);
				ddE(i+n2off,philocs(l)) += (MatP->getdEdkappa())(j) * (ElemP->getddkappadxdphi())(j,i,1,l);
				ddE(i+n3off,philocs(l)) += (MatP->getdEdkappa())(j) * (ElemP->getddkappadxdphi())(j,i,2,l);
			}
		}
	}

	//Symmetries

	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			ddE(philocs(j),i+n1off) = ddE(i+n1off,philocs(j));
			ddE(philocs(j),i+n2off) = ddE(i+n2off,philocs(j));
			ddE(philocs(j),i+n3off) = ddE(i+n3off,philocs(j));
		}
	}

	//dphidphi
	//Second derivative of E(kappa) term
	for (unsigned int i = 0; i < 2; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			for (unsigned int k = 0; k < 3; k++)
			{
				for (unsigned int l = 0; l < 2; l++)
				{
					ddE(philocs(i),philocs(l)) += (ElemP->getdkappadphi())(k,l)*(MatP->getddEddkappa())(j,k)*(ElemP->getdkappadphi())(j,i);
				}
			}
		}
	}
	//Second derivative of kappa term
	for (unsigned int i = 0; i < 2; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			for (unsigned int l = 0; l < 2; l++)
			{
				ddE(philocs(i),philocs(l)) += (MatP->getdEdkappa())(j)*(ElemP->getddkappadphidphi())(j,i,l);
			}
		}
	}

}

void LocalConstructor::setGlobalNodes(Eigen::VectorXi GN)
{
    globalNodes = GN;
}
Eigen::VectorXi LocalConstructor::getGlobalNodes()
{
    return globalNodes;
}


void LocalConstructor::setElementPointer(Element* Ein)
{
	ElemP = Ein;
}

void LocalConstructor::setMaterialPointer(Material* Mat)
{
	MatP = Mat;
}



