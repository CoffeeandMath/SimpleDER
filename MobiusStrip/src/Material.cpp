/*
 * Material.cpp
 *
 *  Created on: Jun 15, 2020
 *      Author: kevin
 */

#include "Material.h"

Material::Material() {
	// TODO Auto-generated constructor stub
	E = 0.0;
	w = 0.1;
	dS = 1.0;
	dSinv = 1.0;
	window = pow(10.0, -5.0);
	dEdkappa.setZero();
	ddEddkappa.setZero();
	nu = 0.3;

	dEdeta.setZero();
	ddEddeta.setZero();
	ddEdkappadeta.setZero();

	kappa0.setZero();
	BMod.setZero();
	Kr.setZero();
	Qr.setZero();

}
void Material::Update(int Order) {
	calcE();
	if (Order > 0) {
		calcdEdkappa();
		calcdEdeta();
	}
	if (Order > 1) {
		calcddEddkappa();
	}
}

void Material::SetElementPointer(Element *Ein) {
	ElemP = Ein;
}
void Material::setBMod(Eigen::VectorXd &Bm) {
	BMod = Bm;
}
double Material::getE() {
	return E;
}
Eigen::Vector3d Material::getdEdkappa() {
	return dEdkappa;
}
Eigen::Matrix3d Material::getddEddkappa() {
	return ddEddkappa;
}
void Material::setkappa0(Eigen::Vector3d &k0) {
	kappa0 = k0;
}

void Material::calcE() {

	E = 0.5 * BMod(0) * pow((ElemP->getkappa())(0), 2.0)
			* pow(1.0 + pow((ElemP->getN2P())->geteta(), 2.0), 2.0)
			* g(etaprime())
			+ 0.5 * BMod(1) * pow((ElemP->getN2P())->geteta(), 2.0)
			+ 0.5 * BMod(2) * pow(etaprime(), 2.0)
			+ 0.5 * BMod(3) * pow(etapprime(), 2.0);
	/*E += BMod(0) * w * (ElemP->getkappa())(0)
	 * (Kr(2, 2)
	 - 2.0 * Kr(0, 2) * (ElemP->getN2P())->geteta() * (1.0 - nu)
	 + Kr(0, 0)
	 * (nu
	 + pow((ElemP->getN2P())->geteta(), 2.0)
	 * (1.0 + nu)));
	 */

	E += BMod(0) * w * (ElemP->getkappa()(0))
			* (Qr(0, 0) * pow((ElemP->getN2P())->geteta(), 2.0)
					- 2.0 * Qr(0, 2) * (ElemP->getN2P())->geteta() + Qr(2, 2));
}
void Material::calcdEdkappa() {
	//for (unsigned int i = 0; i < 3; i++) {
	//	dEdkappa(i) = BMod(i) * ((ElemP->getkappa())(i) - kappa0(i));
	//}
	dEdkappa(0) = BMod(0) * (ElemP->getkappa())(0)
			* pow(1.0 + pow((ElemP->getN2P())->geteta(), 2.0), 2.0)
			* g(etaprime());
	/*
	 dEdkappa(0) += BMod(0) * w
	 * (Kr(2, 2) - 2.0 * Kr(0, 2) * (1.0 - nu)
	 + Kr(0, 0)
	 * (nu
	 + pow((ElemP->getN2P())->geteta(), 2.0)
	 * (1.0 + nu)));
	 * */
	dEdkappa(0) += BMod(0) * w
			* (Qr(0, 0) * pow((ElemP->getN2P())->geteta(), 2.0)
					- 2.0 * Qr(0, 2) * (ElemP->getN2P())->geteta() + Qr(2, 2));

}

void Material::calcdEdeta() {
	dEdeta.setZero();
	dEdeta(1) = 2 * BMod(0) * pow((ElemP->getkappa())(0), 2.0)
			* ((ElemP->getN2P())->geteta())
			* (1.0 + pow(((ElemP->getN2P())->geteta()), 2.0)) * g(etaprime());
	dEdeta += +0.5 * BMod(0) * pow((ElemP->getkappa())(0), 2.0)
			* pow(1.0 + pow(((ElemP->getN2P())->geteta()), 2.0), 2.0)
			* dg(etaprime()) * detaprime();
	dEdeta(1) += BMod(1) * ((ElemP->getN2P())->geteta());
	dEdeta += BMod(2) * etaprime() * detaprime();
	dEdeta += BMod(3) * etapprime() * detapprime();
	/*
	 dEdeta(1) += 2.0 * BMod(0) * w * (ElemP->getkappa())(0)
	 * (Kr(0, 0) * (ElemP->getN2P())->geteta() * (1.0 + nu)
	 - Kr(0, 2) * (1.0 - nu));
	 */
	dEdeta(1) += 2.0 * BMod(0) * w * (ElemP->getkappa())(0)
			* (((ElemP->getN2P())->geteta()) * Qr(0, 0) - Qr(0, 2));

}

void Material::calcddEddkappa() {
	//ddEddkappa.setZero();
	double gval = g(etaprime());
	Eigen::Vector3d dgdeta = dg(etaprime()) * detaprime();
	Eigen::Matrix3d ddgddeta = dg(etaprime()) * ddetaprime()
			+ ddg(etaprime()) * detaprime() * (detaprime().transpose());

	double dfdkappa1 = BMod(0) * (ElemP->getkappa())(0)
			* pow(1.0 + pow(((ElemP->getN2P())->geteta()), 2.0), 2.0);
	Eigen::Vector3d dfdeta;
	dfdeta << 0.0, 2.0 * BMod(0) * pow((ElemP->getkappa())(0), 2.0)
			* ((ElemP->getN2P())->geteta())
			* (1.0 + pow(((ElemP->getN2P())->geteta()), 2.0)), 0.0;
	double ddfddkappa1 = BMod(0)
			* pow(1.0 + pow(((ElemP->getN2P())->geteta()), 2.0), 2.0);

	Eigen::Vector3d ddfdkappa1deta;
	ddfdkappa1deta << 0.0, BMod(0) * 4.0 * (ElemP->getkappa())(0)
			* ((ElemP->getN2P())->geteta())
			* (1.0 + pow(((ElemP->getN2P())->geteta()), 2.0)), 0.0;

	Eigen::Matrix3d ddfddeta;
	ddfddeta << 0.0, 0.0, 0.0, 0.0, BMod(0) * 2.0
			* pow((ElemP->getkappa())(0), 2.0)
			* (1.0 + 3.0 * pow(((ElemP->getN2P())->geteta()), 2.0)), 0.0, 0.0, 0.0, 0.0;

	ddEddkappa << ddfddkappa1 * gval, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	Eigen::Vector3d vtmp;
	vtmp = gval * ddfdkappa1deta + dfdkappa1 * dgdeta;
	ddEdkappadeta << vtmp(0), vtmp(1), vtmp(2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	ddEddeta = ddfddeta * gval + dfdeta * dgdeta.transpose()
			+ dgdeta * dfdeta.transpose()
			+ 0.5 * BMod(0) * pow((ElemP->getkappa())(0), 2.0)
					* pow(1.0 + pow((ElemP->getN2P())->geteta(), 2.0), 2.0)
					* ddgddeta + BMod(2) * detaprime() * detaprime().transpose()
			+ BMod(3) * detapprime() * detapprime().transpose();
	ddEddeta(1, 1) += BMod(1);

	/*
	 ddEdkappadeta(0, 1) += 2.0 * w * BMod(0)
	 * (Kr(0, 0) * ((ElemP->getN2P())->geteta()) * (1.0 + nu)
	 - (1.0 - nu) * Kr(0, 2));
	 ddEddeta(1, 1) += 2.0 * BMod(0) * w * Kr(0, 0) * (1.0 + nu)
	 * (ElemP->getkappa())(0);
	 * */
	ddEdkappadeta(0, 1) += 2.0 * BMod(0) * w
			* (((ElemP->getN2P())->geteta()) * Qr(0, 0) - Qr(0, 2));
	ddEddeta(1, 1) += 2.0 * BMod(0) * w * Qr(0, 0) * (ElemP->getkappa())(0);
}

double Material::g(double ep) {
	double gval = 0.0;

	
	if (abs(ep) < window) {
		gval =  w + pow(w, 3.0) * pow(ep, 2.0) / 12.0
				+ pow(w, 5.0) * pow(ep, 4.0) / 80.0;
	} else {
		gval =  log((2.0 + ep * w) / (2.0 - ep * w)) / ep;
		/*
		 
		 return w + pow(w, 3.0) * pow(ep, 2.0)
		 + pow(w, 5.0) * pow(ep, 4.0) / 80.0
		 + pow(w, 7.0) * pow(ep, 6.0) / 448.0
		 + pow(w, 9.0) * pow(ep, 8.0) / 2304.0;
		 */
	}



	if (isnan(gval)){
		//std::cout << "ep too large" << std::endl;
	}
	return gval;


}
double Material::dg(double ep) {
	if (abs(ep) < window) {
		return pow(w, 3.0) * ep / 6.0 + pow(w, 5.0) * pow(ep, 3.0) / 20.0;
	} else {
		return ((4.0 * w * ep) / (4.0 - pow(w * ep, 2.0))
				- log(-1.0 - 4.0 / (w * ep - 2.0))) / pow(ep, 2.0);
		/*
		 return pow(w, 3.0) * ep / 6.0 + pow(w, 5.0) * pow(ep, 3.0) / 20.0
		 + 3.0 * pow(w, 7.0) * pow(ep, 5.0) / 244.0
		 + pow(w, 9.0) * pow(ep, 7.0) / 288.0;
		 */
	}

}
double Material::ddg(double ep) {
	if (abs(ep) < window) {
		return pow(w, 3.0) / 6.0 + 3.0 * pow(w, 5.0) * pow(ep, 2.0) / 20.0
				+ 15.0 * pow(w, 7.0) * pow(ep, 4.0) / 224.0;
	} else {
		return 2.0
				* (8.0 * ep * w * (pow(ep, 2.0) * pow(w, 2.0) - 2.0)
						/ pow(pow(ep * w, 2.0) - 4.0, 2.0)
						+ log(-1.0 - 4.0 / (ep * w - 2.0))) / pow(ep, 3.0);
		/*
		 return pow(w, 3.0) / 6.0 + 3.0 * pow(w, 5.0) * pow(ep, 2.0) / 20.0
		 + 15.0 * pow(w, 7.0) * pow(ep, 4.0) / 244.0
		 + 7.0 * pow(w, 9.0) * pow(ep, 6.0) / 288.0
		 + 45.0 * pow(w, 11.0) * pow(ep, 8.0) / 5632.0;
		 */

	}
}
void Material::setdS(double dStemp) {
	dS = dStemp;
	dSinv = 1.0 / dStemp;
}
double Material::etaprime() {
	//return 0.5 * dSinv
	//		* ((ElemP->getN3P())->geteta() - (ElemP->getN1P())->geteta());
	double etapval = dSinv*((ElemP->getN3P())->geteta() - (ElemP->getN2P())->geteta());
	/*
	double evunsc = 1.9999;
	if (abs(w*etapval) >= 2.0) {
		if (etapval > 0.0){
			etapval = evunsc/w;
		} else {
			etapval = -1.0*evunsc/w;
		}
	}
	*/
	return etapval;
}

double Material::etapprime() {
	return pow(dSinv, 2.0)
			* ((ElemP->getN3P())->geteta() - 2.0 * (ElemP->getN2P())->geteta()
					+ (ElemP->getN1P())->geteta());
}

Eigen::Vector3d Material::detapprime() {

	double sc = pow(dSinv, 2.0);
	Eigen::Vector3d tmp = { sc, -2.0 * sc, sc };
	return tmp;
}
Eigen::Vector3d Material::detaprime() {
	//Eigen::Vector3d tmp = { -0.5 * dSinv, 0.0, 0.5 * dSinv };

	Eigen::Vector3d tmp = {0.0,-1.0*dSinv,dSinv};
	return tmp;
}
Eigen::Matrix3d Material::ddetaprime() {
	Eigen::Matrix3d tmp;
	tmp.setZero();
	return tmp;
}
Eigen::Vector3d Material::getdEdeta() {
	return dEdeta;
}
Eigen::Matrix3d Material::getddEdkappadeta() {
	return ddEdkappadeta;
}
Eigen::Matrix3d Material::getddEddeta() {
	return ddEddeta;
}
void Material::setw(double wt) {
	w = wt;
}

void Material::setKr(Eigen::Matrix3d &Krtemp) {
	Kr = Krtemp;
	Qr = Kr + nu * cof(Kr);
}
void Material::setnu(double n) {
	nu = n;
}
Eigen::Matrix3d Material::getKr() {
	return Kr;
}

Eigen::Matrix3d Material::cof(Eigen::Matrix3d &F) {
	Eigen::Matrix3d cofF;
	cofF.setZero();
	cofF(0, 0) = F(2, 2);
	cofF(0, 2) = -F(0, 2);
	cofF(2, 0) = -F(2, 0);
	cofF(2, 2) = F(0, 0);
	return cofF;
}

double Material::getetaprime(){
	double etat = etapprime();
	return etat;

}
Eigen::Vector3d Material::getdetaprime(){
	Eigen::Vector3d detat;
	detat = detaprime();
	return detat;

}
Eigen::Matrix3d Material::getddetaprime(){
	Eigen::Matrix3d ddetat;
	ddetat = ddetaprime();
	return ddetat;

}
double Material::getw(){
	return w;
}