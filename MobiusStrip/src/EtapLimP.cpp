/*
 * EtapLimP.cpp
 *
 *  Created on: November 17, 2020
 *      Author: KevinKorner
 */

#include "EtapLimP.h"

EtapLimP::EtapLimP() {
	MatP = NULL;
	dC.resize(3);
	dC.setZero();
	ddC.resize(3, 3);
	ddC.setZero();
	Bloc = 1.99;
	
}

void EtapLimP::Update(int Order) {
	calcC();
	if (Order > 0) {
		calcdC();
	}
	if (Order > 1) {
		calcddC();
	}
}

void EtapLimP::calcC() {
	C = (MatP->getetaprime())*(MatP->getw()) + Bloc;
}
void EtapLimP::calcdC() {
	Eigen::Vector3d detatemp = (MatP->getdetaprime());
	double wtemp = (MatP->getw());
	for (unsigned int i = 0; i < 3; i++){
		dC[i] =  detatemp[i]* wtemp;
	}

	

}
void EtapLimP::calcddC() {
	Eigen::Matrix3d ddetatemp = (MatP->getddetaprime());
	double wtemp = (MatP->getw());
	for (unsigned int i = 0; i < 3; i++){
		for (unsigned int j = 0; j < 3; j++){
			ddC(i,j) =  ddetatemp(i,j)* wtemp;
		}
	}
}
void EtapLimP::setMaterialPointer(Material *np) {
	MatP = np;
}

void EtapLimP::setBarrierLoc(double bl) {
	Bloc = bl;
}