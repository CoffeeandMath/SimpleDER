/*
 * EtapLimN.cpp
 *
 *  Created on: November 17, 2020
 *      Author: KevinKorner
 */

#include "EtapLimN.h"

EtapLimN::EtapLimN() {
	MatP = NULL;
	dC.resize(3);
	dC.setZero();
	ddC.resize(3, 3);
	ddC.setZero();
	Bloc = 1.99;
	
}

void EtapLimN::Update(int Order) {
	calcC();
	if (Order > 0) {
		calcdC();
	}
	if (Order > 1) {
		calcddC();
	}
}

void EtapLimN::calcC() {
	C = Bloc - (MatP->getetaprime())*(MatP->getw());
}
void EtapLimN::calcdC() {
	Eigen::Vector3d detatemp = (MatP->getdetaprime());
	double wtemp = (MatP->getw());
	for (unsigned int i = 0; i < 3; i++){
		dC[i] =  -1.0*detatemp[i]* wtemp;
	}

	

}
void EtapLimN::calcddC() {
	Eigen::Matrix3d ddetatemp = (MatP->getddetaprime());
	double wtemp = (MatP->getw());
	for (unsigned int i = 0; i < 3; i++){
		for (unsigned int j = 0; j < 3; j++){
			ddC(i,j) =  -1.0*ddetatemp(i,j)* wtemp;
		}
	}


	

}
void EtapLimN::setMaterialPointer(Material *np) {
	MatP = np;
}

void EtapLimN::setBarrierLoc(double bl) {
	Bloc = bl;
}