/*
 * RibbonConst.cpp
 *
 *  Created on: July 1, 2020
 *      Author: KevinKorner
 */

#include "RibbonConst.h"

RibbonConst::RibbonConst() {
	// TODO Auto-generated constructor stub
	dS = 1.0;
	dSinv = 1.0;
	dC.resize(14);
	dC.setZero();
	ddC.resize(14, 14);
	ddC.setZero();

}
void RibbonConst::Update(int Order) {
	calcC();
	if (Order > 0) {
		calcdC();
	}
	if (Order > 1) {
		calcddC();
	}
}
void RibbonConst::setElementPointer(Element *ep) {
	ElemP = ep;
}
void RibbonConst::setdS(double dst) {
	dS = dst;
	dSinv = 1.0 / dst;
}
void RibbonConst::calcC() {
	C = dSinv
			* ((ElemP->getkappa())(2)
					- ((ElemP->getN2P())->geteta()) * ((ElemP->getkappa())(0)));
}
void RibbonConst::calcdC() {

	Eigen::VectorXd dk1dq(14);
	dk1dq << (ElemP->getdkappadx())(0, 0, 0), (ElemP->getdkappadx())(0, 1, 0), (ElemP->getdkappadx())(
			0, 2, 0), 0.0, (ElemP->getdkappadphi())(0, 0), (ElemP->getdkappadx())(
			0, 0, 1), (ElemP->getdkappadx())(0, 1, 1), (ElemP->getdkappadx())(0,
			2, 1), 0.0, (ElemP->getdkappadphi())(0, 1), (ElemP->getdkappadx())(
			0, 0, 2), (ElemP->getdkappadx())(0, 1, 2), (ElemP->getdkappadx())(0,
			2, 2), 0.0;

	Eigen::VectorXd dk3dq(14);
	dk3dq << (ElemP->getdkappadx())(2, 0, 0), (ElemP->getdkappadx())(2, 1, 0), (ElemP->getdkappadx())(
			2, 2, 0), 0.0, (ElemP->getdkappadphi())(2, 0), (ElemP->getdkappadx())(
			2, 0, 1), (ElemP->getdkappadx())(2, 1, 1), (ElemP->getdkappadx())(2,
			2, 1), 0.0, (ElemP->getdkappadphi())(2, 1), (ElemP->getdkappadx())(
			2, 0, 2), (ElemP->getdkappadx())(2, 1, 2), (ElemP->getdkappadx())(2,
			2, 2), 0.0;

	dC = dSinv * (dk3dq - ((ElemP->getN2P())->geteta()) * dk1dq);
	dC(8) = -dSinv * (ElemP->getkappa())(0);
}
void RibbonConst::calcddC() {
	int n1off = 0;
	int n2off = 5;
	int n3off = 10;
	Eigen::Vector2i philocs;
	philocs << 4, 9;
	Eigen::Vector3i etalocs;
	etalocs << 3, 8, 13;

	Eigen::MatrixXd ddk1ddq(14, 14);
	ddk1ddq.setZero();

	Eigen::MatrixXd ddk3ddq(14, 14);
	ddk3ddq.setZero();

	for (unsigned int i = 0; i < 3; i++) {
		for (unsigned int j = 0; j < 3; j++) {
			ddk1ddq(i + n1off, j + n1off) = dSinv
					* (ElemP->getddkappadxdx())(0, i, j, 0, 0);
			ddk1ddq(i + n1off, j + n2off) = dSinv
					* (ElemP->getddkappadxdx())(0, i, j, 0, 1);
			ddk1ddq(i + n1off, j + n3off) = dSinv
					* (ElemP->getddkappadxdx())(0, i, j, 0, 2);
			ddk1ddq(i + n2off, j + n1off) = dSinv
					* (ElemP->getddkappadxdx())(0, i, j, 1, 0);
			ddk1ddq(i + n2off, j + n2off) = dSinv
					* (ElemP->getddkappadxdx())(0, i, j, 1, 1);
			ddk1ddq(i + n2off, j + n3off) = dSinv
					* (ElemP->getddkappadxdx())(0, i, j, 1, 2);
			ddk1ddq(i + n3off, j + n1off) = dSinv
					* (ElemP->getddkappadxdx())(0, i, j, 2, 0);
			ddk1ddq(i + n3off, j + n2off) = dSinv
					* (ElemP->getddkappadxdx())(0, i, j, 2, 1);
			ddk1ddq(i + n3off, j + n3off) = dSinv
					* (ElemP->getddkappadxdx())(0, i, j, 2, 2);

			ddk3ddq(i + n1off, j + n1off) = dSinv
					* (ElemP->getddkappadxdx())(2, i, j, 0, 0);
			ddk3ddq(i + n1off, j + n2off) = dSinv
					* (ElemP->getddkappadxdx())(2, i, j, 0, 1);
			ddk3ddq(i + n1off, j + n3off) = dSinv
					* (ElemP->getddkappadxdx())(2, i, j, 0, 2);
			ddk3ddq(i + n2off, j + n1off) = dSinv
					* (ElemP->getddkappadxdx())(2, i, j, 1, 0);
			ddk3ddq(i + n2off, j + n2off) = dSinv
					* (ElemP->getddkappadxdx())(2, i, j, 1, 1);
			ddk3ddq(i + n2off, j + n3off) = dSinv
					* (ElemP->getddkappadxdx())(2, i, j, 1, 2);
			ddk3ddq(i + n3off, j + n1off) = dSinv
					* (ElemP->getddkappadxdx())(2, i, j, 2, 0);
			ddk3ddq(i + n3off, j + n2off) = dSinv
					* (ElemP->getddkappadxdx())(2, i, j, 2, 1);
			ddk3ddq(i + n3off, j + n3off) = dSinv
					* (ElemP->getddkappadxdx())(2, i, j, 2, 2);
		}
	}

	for (unsigned int i = 0; i < 3; i++) {
		for (unsigned int j = 0; j < 2; j++) {
			ddk1ddq(i + n1off, philocs(j)) = dSinv
					* (ElemP->getddkappadxdphi())(0, i, 0, j);
			ddk1ddq(i + n2off, philocs(j)) = dSinv
					* (ElemP->getddkappadxdphi())(0, i, 1, j);
			ddk1ddq(i + n3off, philocs(j)) = dSinv
					* (ElemP->getddkappadxdphi())(0, i, 2, j);

			ddk1ddq(philocs(j), i + n1off) = dSinv
					* (ElemP->getddkappadxdphi())(0, i, 0, j);
			ddk1ddq(philocs(j), i + n2off) = dSinv
					* (ElemP->getddkappadxdphi())(0, i, 1, j);
			ddk1ddq(philocs(j), i + n3off) = dSinv
					* (ElemP->getddkappadxdphi())(0, i, 2, j);

			ddk3ddq(i + n1off, philocs(j)) = dSinv
					* (ElemP->getddkappadxdphi())(2, i, 0, j);
			ddk3ddq(i + n2off, philocs(j)) = dSinv
					* (ElemP->getddkappadxdphi())(2, i, 1, j);
			ddk3ddq(i + n3off, philocs(j)) = dSinv
					* (ElemP->getddkappadxdphi())(2, i, 2, j);

			ddk3ddq(philocs(j), i + n1off) = dSinv
					* (ElemP->getddkappadxdphi())(2, i, 0, j);
			ddk3ddq(philocs(j), i + n2off) = dSinv
					* (ElemP->getddkappadxdphi())(2, i, 1, j);
			ddk3ddq(philocs(j), i + n3off) = dSinv
					* (ElemP->getddkappadxdphi())(2, i, 2, j);
		}
	}

	for (unsigned int i = 0; i < 2; i++) {
		for (unsigned int j = 0; j < 2; j++) {
			ddk1ddq(philocs(i), philocs(j)) = dSinv
					* (ElemP->getddkappadphidphi())(0, i, j);

			ddk3ddq(philocs(i), philocs(j)) = dSinv
					* (ElemP->getddkappadphidphi())(2, i, j);
		}
	}

	Eigen::VectorXd dk1dq(14);
	dk1dq << dSinv * (ElemP->getdkappadx())(0, 0, 0), dSinv
			* (ElemP->getdkappadx())(0, 1, 0), dSinv
			* (ElemP->getdkappadx())(0, 2, 0), 0.0, dSinv
			* (ElemP->getdkappadphi())(0, 0), dSinv
			* (ElemP->getdkappadx())(0, 0, 1), dSinv
			* (ElemP->getdkappadx())(0, 1, 1), dSinv
			* (ElemP->getdkappadx())(0, 2, 1), 0.0, dSinv
			* (ElemP->getdkappadphi())(0, 1), dSinv
			* (ElemP->getdkappadx())(0, 0, 2), dSinv
			* (ElemP->getdkappadx())(0, 1, 2), dSinv
			* (ElemP->getdkappadx())(0, 2, 2), 0.0;

	ddC = ddk3ddq - ((ElemP->getN2P())->geteta()) * ddk1ddq;
	for (unsigned int i = 0; i < 14; i++) {
		ddC(i, 8) = -dk1dq(i);
		ddC(8, i) = -dk1dq(i);
	}

}

