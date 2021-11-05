/*
 * Element.cpp
 *
 *  Created on: Jun 12, 2020
 *      Author: KevinKorner
 */

#include "Element.h"

Element::Element() {
	// TODO Auto-generated constructor stub
	q.Set(0.0,0.0,0.0,0.0);
	kappa.setZero();

	dkappadphi.setZero();
	dkappadx = Eigen::Tensor<double, 3> (3,3,3);
	dkappadx.setZero();
	ddkappadphidphi = Eigen::Tensor<double, 3> (3,2,2);
	ddkappadphidphi.setZero();
	ddkappadxdphi = Eigen::Tensor<double, 4> (3,3,3,2);
	ddkappadxdphi.setZero();
	ddkappadxdx = Eigen::Tensor<double, 5> (3,3,3,3,3);
	ddkappadxdx.setZero();
	qe3.Set(0.0,0.0,0.0,1.0);


}

void Element::setFramesPointers(Frame* Frame1,Frame* Frame2)
{
	Frame1P = Frame1;
	Frame2P = Frame2;
}

void Element::setNodesPointers(Node* N1,Node* N2,Node* N3)
{
	Node1P = N1;
	Node2P = N2;
	Node3P = N3;
}

void Element::Update(int Order)
{
	q = (Frame1P->getd()).conj() * (Frame2P->getd());
	kappa = 2.0*q.Getv();

	if (Order > 0)
	{

		calcdkappadx();
		calcdkappadphi();
	}
	if (Order > 1)
	{
		calcddkappadphidphi();
		calcddkappadxdphi();
		calcddkappadxdx();
	}
}


void Element::calcdkappadx()
{
	for (unsigned int i = 0; i < 3; i++)
	{
		Quat<double> Bimquat;
		Bimquat.Set(0.0,(Frame1P->getB())(0,i),(Frame1P->getB())(1,i),(Frame1P->getB())(2,i));
		//Bimquat.Set(0.0, Bim(0,i),Bim(1,i),Bim(2,i));
		Quat<double> temp = 2.0*(Frame1P->getd()).conj() * Bimquat * (Frame2P->getd());
		dkappadx(0,i,0) = temp.Getx();
		dkappadx(1,i,0) = temp.Gety();
		dkappadx(2,i,0) = temp.Getz();
	}
	for (unsigned int i = 0; i < 3; i++)
		{
			Quat<double> Biquat;
			Biquat.Set(0.0,(Frame2P->getB())(0,i),(Frame2P->getB())(1,i),(Frame2P->getB())(2,i));
			//Biquat.Set(0.0, Bi(0,i),Bi(1,i),Bi(2,i));
			Quat<double> temp = 2.0*(Frame1P->getd()).conj() * Biquat * (Frame2P->getd());
			dkappadx(0,i,2) = temp.Getx();
			dkappadx(1,i,2) = temp.Gety();
			dkappadx(2,i,2) = temp.Getz();
		}
	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			dkappadx(i,j,1) = -dkappadx(i,j,2) - dkappadx(i,j,0);
		}
	}
}

void Element::calcdkappadphi()
{
	Quat<double> temp;
	temp = -1.0*qe3*q;
	dkappadphi(0,0) = temp.Getx();
	dkappadphi(1,0) = temp.Gety();
	dkappadphi(2,0) = temp.Getz();

	temp = q*qe3;
	dkappadphi(0,1) = temp.Getx();
	dkappadphi(1,1) = temp.Gety();
	dkappadphi(2,1) = temp.Getz();



}

void Element::calcddkappadphidphi()
{
	Quat<double> qtemp = qe3*q*qe3;

	for (unsigned int i = 0; i < 3; i++)
	{
		ddkappadphidphi(i,0,0) = -0.5*q.Getind(i+1);
		ddkappadphidphi(i,0,1) = -0.5*qtemp.Getind(i+1);
		ddkappadphidphi(i,1,0) = ddkappadphidphi(i,0,1);
		ddkappadphidphi(i,1,1) = ddkappadphidphi(i,0,0);
	}
}

void Element::calcddkappadxdphi()
{
	for (unsigned int i = 0; i < 3; i++)
		{
			Quat<double> Bimquat;
			Quat<double> Biquat;

			Bimquat.Set(0.0,(Frame1P->getB())(0,i),(Frame1P->getB())(1,i),(Frame1P->getB())(2,i));
			Biquat.Set(0.0,(Frame2P->getB())(0,i),(Frame2P->getB())(1,i),(Frame2P->getB())(2,i));
			//Biquat.Set(0.0, Bi(0,i),Bi(1,i),Bi(2,i));
			Quat<double> dqndphinmdxnm = -1.0*qe3*(Frame1P->getd()).conj()*Bimquat*(Frame2P->getd());
			Quat<double> dqndphinmdxnp = -1.0*qe3*(Frame1P->getd()).conj()*Biquat*(Frame2P->getd());
			Quat<double> dqndphindxnm = (Frame1P->getd()).conj()*Bimquat*(Frame2P->getd())*qe3;
			Quat<double> dqndphindxnp = (Frame1P->getd()).conj()*Biquat*(Frame2P->getd())*qe3;
			for (unsigned int j = 0; j < 3; j++)
			{
				ddkappadxdphi(j,i,0,0) = dqndphinmdxnm.Getind(j+1);
				ddkappadxdphi(j,i,2,0) = dqndphinmdxnp.Getind(j+1);
				ddkappadxdphi(j,i,1,0) = -(1.0)*ddkappadxdphi(j,i,2,0) - ddkappadxdphi(j,i,0,0);


				ddkappadxdphi(j,i,0,1) = dqndphindxnm.Getind(j+1);
				ddkappadxdphi(j,i,2,1) = dqndphindxnp.Getind(j+1);
				ddkappadxdphi(j,i,1,1) = -(1.0)*ddkappadxdphi(j,i,0,1) - ddkappadxdphi(j,i,2,1);
			}

		}

}

void Element::calcddkappadxdx()
{
	Quat<double> cdim;
	cdim = (Frame1P->getd()).conj();

	for (unsigned int n = 0; n < 3; n++)
	{
		Quat<double> Bimquatn;
		Quat<double> Biquatn;
		Bimquatn.Set(0.0,(Frame1P->getB())(0,n),(Frame1P->getB())(1,n),(Frame1P->getB())(2,n));
		Biquatn.Set(0.0,(Frame2P->getB())(0,n),(Frame2P->getB())(1,n),(Frame2P->getB())(2,n));
		for (unsigned int q = 0; q < 3; q++)
		{
			Quat<double> Bimquatq;
			Quat<double> Biquatq;
			Bimquatq.Set(0.0,(Frame1P->getB())(0,q),(Frame1P->getB())(1,q),(Frame1P->getB())(2,q));
			Biquatq.Set(0.0,(Frame2P->getB())(0,q),(Frame2P->getB())(1,q),(Frame2P->getB())(2,q));

			Quat<double> DBimnq;
			DBimnq.Set(0.0,(Frame1P->getDB())(0,n,q),(Frame1P->getDB())(1,n,q),(Frame1P->getDB())(2,n,q));

			Quat<double> DBinq;
			DBinq.Set(0.0,(Frame2P->getDB())(0,n,q),(Frame2P->getDB())(1,n,q),(Frame2P->getDB())(2,n,q));

			Quat<double> temp1;
			Quat<double> temp2;
			Quat<double> temp3;

			temp1 = cdim*(Bimquatq*Bimquatn -DBimnq)*(Frame2P->getd());
			temp2 = cdim*Bimquatn*Biquatq*(Frame2P->getd());
			temp3 = cdim*(Biquatn*Biquatq + DBinq)*(Frame2P->getd());

			for (unsigned int k = 0; k<3; k++)
			{
				ddkappadxdx(k,n,q,0,0) = 2.0*temp1.Getind(k+1);
				ddkappadxdx(k,n,q,0,2) = 2.0*temp2.Getind(k+1);
				ddkappadxdx(k,n,q,2,2) = 2.0*temp3.Getind(k+1);
			}

		}
	}

	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j<3; j++)
		{
			for (unsigned int k = 0; k < 3; k++)
			{
				ddkappadxdx(i,j,k,2,0) = ddkappadxdx(i,k,j,0,2);
			}
		}
	}

	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j<3; j++)
		{
			for (unsigned int k = 0; k < 3; k++)
			{
				ddkappadxdx(i,j,k,0,1) = -1.0*ddkappadxdx(i,j,k,0,2) -ddkappadxdx(i,j,k,0,0);
			}
		}
	}

	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j<3; j++)
		{
			for (unsigned int k = 0; k < 3; k++)
			{
				ddkappadxdx(i,j,k,1,0) = ddkappadxdx(i,k,j,0,1);
			}
		}
	}

	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j<3; j++)
		{
			for (unsigned int k = 0; k < 3; k++)
			{
				ddkappadxdx(i,j,k,1,2) = -1.0*ddkappadxdx(i,j,k,0,2) -ddkappadxdx(i,j,k,2,2);
			}
		}
	}
	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j<3; j++)
		{
			for (unsigned int k = 0; k < 3; k++)
			{
				ddkappadxdx(i,j,k,2,1) = ddkappadxdx(i,k,j,1,2);
			}
		}
	}

	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j<3; j++)
		{
			for (unsigned int k = 0; k < 3; k++)
			{
				ddkappadxdx(i,j,k,1,1) = -1.0*ddkappadxdx(i,j,k,1,2) -ddkappadxdx(i,j,k,1,0);
			}
		}
	}

}

void Element::setGlobalNodes(Eigen::VectorXi GN)
{
    globalNodes = GN;
}
Eigen::VectorXi Element::getGlobalNodes()
{
    return globalNodes;
}

Quat<double> Element::getq()
{
	return q;
}
Eigen::Vector3d Element::getkappa()
{
	return kappa;
}
Eigen::Matrix<double,3,2> Element::getdkappadphi()
{
	return dkappadphi;
}
Eigen::Tensor<double,3> Element::getdkappadx()
{
	return dkappadx;
}
Eigen::Tensor<double,3> Element::getddkappadphidphi()
{
	return ddkappadphidphi;
}
Eigen::Tensor<double,4> Element::getddkappadxdphi()
{
	return ddkappadxdphi;
}
Eigen::Tensor<double,5> Element::getddkappadxdx()
{
	return ddkappadxdx;
}
