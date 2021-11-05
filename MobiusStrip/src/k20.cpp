/*
 * k20.cpp
 *
 *  Created on: Jun 30, 2020
 *      Author: KevinKorner
 */

#include "k20.h"

k20::k20()
{
	// TODO Auto-generated constructor stub
	dS = 1.0;
	dSinv = 1.0;
	dC.resize(14);
	dC.setZero();
	ddC.resize(14, 14);
	ddC.setZero();

}
void k20::Update(int Order)
{
	calcC();
	if (Order > 0)
	{
		calcdC();
	}
	if (Order > 1)
	{
		calcddC();
	}
}
void k20::setElementPointer(Element *ep)
{
	ElemP = ep;
}
void k20::setdS(double dst)
{
	dS = dst;
	dSinv = 1.0 / dst;
}
void k20::calcC()
{
	C = dSinv * (ElemP->getkappa())(1);
}
void k20::calcdC()
{

	dC << dSinv * (ElemP->getdkappadx())(1, 0, 0), dSinv
			* (ElemP->getdkappadx())(1, 1, 0), dSinv
			* (ElemP->getdkappadx())(1, 2, 0), 0.0, dSinv
			* (ElemP->getdkappadphi())(1, 0), dSinv
			* (ElemP->getdkappadx())(1, 0, 1), dSinv
			* (ElemP->getdkappadx())(1, 1, 1), dSinv
			* (ElemP->getdkappadx())(1, 2, 1), 0.0, dSinv
			* (ElemP->getdkappadphi())(1, 1), dSinv
			* (ElemP->getdkappadx())(1, 0, 2), dSinv
			* (ElemP->getdkappadx())(1, 1, 2), dSinv
			* (ElemP->getdkappadx())(1, 2, 2), 0.0;
}
void k20::calcddC()
{
	int n1off = 0;
	int n2off = 5;
	int n3off = 10;
	Eigen::Vector2i philocs;
	philocs << 4, 9;
	Eigen::Vector3i etalocs;
	etalocs << 3, 8, 13;

	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			ddC(i + n1off, j + n1off) = dSinv
					* (ElemP->getddkappadxdx())(1, i, j, 0, 0);
			ddC(i + n1off, j + n2off) = dSinv
					* (ElemP->getddkappadxdx())(1, i, j, 0, 1);
			ddC(i + n1off, j + n3off) = dSinv
					* (ElemP->getddkappadxdx())(1, i, j, 0, 2);
			ddC(i + n2off, j + n1off) = dSinv
					* (ElemP->getddkappadxdx())(1, i, j, 1, 0);
			ddC(i + n2off, j + n2off) = dSinv
					* (ElemP->getddkappadxdx())(1, i, j, 1, 1);
			ddC(i + n2off, j + n3off) = dSinv
					* (ElemP->getddkappadxdx())(1, i, j, 1, 2);
			ddC(i + n3off, j + n1off) = dSinv
					* (ElemP->getddkappadxdx())(1, i, j, 2, 0);
			ddC(i + n3off, j + n2off) = dSinv
					* (ElemP->getddkappadxdx())(1, i, j, 2, 1);
			ddC(i + n3off, j + n3off) = dSinv
					* (ElemP->getddkappadxdx())(1, i, j, 2, 2);
		}
	}

	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			ddC(i + n1off, philocs(j)) = dSinv
					* (ElemP->getddkappadxdphi())(1, i, 0, j);
			ddC(i + n2off, philocs(j)) = dSinv
					* (ElemP->getddkappadxdphi())(1, i, 1, j);
			ddC(i + n3off, philocs(j)) = dSinv
					* (ElemP->getddkappadxdphi())(1, i, 2, j);

			ddC(philocs(j), i + n1off) = dSinv
					* (ElemP->getddkappadxdphi())(1, i, 0, j);
			ddC(philocs(j), i + n2off) = dSinv
					* (ElemP->getddkappadxdphi())(1, i, 1, j);
			ddC(philocs(j), i + n3off) = dSinv
					* (ElemP->getddkappadxdphi())(1, i, 2, j);

		}
	}

	for (unsigned int i = 0; i < 2; i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			ddC(philocs(i), philocs(j)) = dSinv
					* (ElemP->getddkappadphidphi())(1, i, j);
		}
	}
}

