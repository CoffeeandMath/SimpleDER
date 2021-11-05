/*
 * Inextensibility.cpp
 *
 *  Created on: Jun 19, 2020
 *      Author: KevinKorner
 */

#include "Inextensibility.h"

Inextensibility::Inextensibility() {
	// TODO Auto-generated constructor stub
	FrameP = NULL;
	dS = 1.0;
	dSinv = 1.0;
	dC.resize(7);
	dC.setZero();
    ddC.resize(7,7);
    ddC.setZero();
}
void Inextensibility::Update(int Order)
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
void Inextensibility::calcC()
{
	C = (FrameP->geteps())/dS - 1.0;
}
void Inextensibility::calcdC()
{

	for (unsigned int i = 0; i < 3; i++)
	{
		dC[i] = dSinv*(FrameP->getdepsdx())(i,0);
		dC[i+4] = dSinv*(FrameP->getdepsdx())(i,1);
	}

}
void Inextensibility::calcddC()
{
	int n1off = 0;
	int n2off = 4;
	for (unsigned int i = 0; i < 3; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			ddC(i+n1off,j+n1off) = dSinv*(FrameP->getddepsdxdx())(i,j,0,0);
			ddC(i+n1off,j+n2off) = dSinv*(FrameP->getddepsdxdx())(i,j,0,1);
			ddC(i+n2off,j+n1off) = dSinv*(FrameP->getddepsdxdx())(i,j,1,0);
			ddC(i+n2off,j+n2off) = dSinv*(FrameP->getddepsdxdx())(i,j,1,1);
		}
	}

}
void Inextensibility::setFramePointer(Frame* fp)
{
	FrameP = fp;
}
void Inextensibility::setdS(double d)
{
	dS = d;
	dSinv = 1.0/dS;
}
