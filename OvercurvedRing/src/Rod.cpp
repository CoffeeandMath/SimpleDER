/*
 * Rod.cpp
 *
 *  Created on: Jun 16, 2020
 *      Author: kevin
 */

#include "Rod.h"

Rod::Rod() {
	// TODO Auto-generated constructor stub
	Nn = 0;
	Parallelized = false;
	E = 0.0;
	Ebend = 0.0;
	Egrav = 0.0;
	Eext = 0.0;
	gravdir << 0.0, -1.0, -0.0;
	NvecP = NULL;
	FvecP = NULL;
	EvecP = NULL;
	MvecP = NULL;
	LCvecP = NULL;
	rho = 0.0;
	Ndof = 0;
	PConstP = NULL;
	Nc = 0;

}

void Rod::Update(int Order) {
	if (Parallelized) {
#pragma omp parallel
		{
#pragma omp for
			for (unsigned int i = 0; i < (*FvecP).size(); i++) {
				((*FvecP)[i]).Update(Order);
			}

#pragma omp for
			for (unsigned int i = 0; i < (*EvecP).size(); i++) {
				((*EvecP)[i]).Update(Order);

				((*MvecP)[i]).Update(Order);

				((*LCvecP)[i]).Update(Order);
			}

			if (PConstP != NULL) {
#pragma omp for
				for (unsigned int i = 0; i < PConstP->size(); i++) {
					((*PConstP)[i])->Update(Order);
				}
			}

		}

	} else {

		for (unsigned int i = 0; i < (*FvecP).size(); i++) {
			((*FvecP)[i]).Update(Order);
		}

		for (unsigned int i = 0; i < (*EvecP).size(); i++) {
			((*EvecP)[i]).Update(Order);
		}

		for (unsigned int i = 0; i < (*MvecP).size(); i++) {
			((*MvecP)[i]).Update(Order);
		}

		for (unsigned int i = 0; i < (*LCvecP).size(); i++) {
			((*LCvecP)[i]).Update(Order);
		}

		if (PConstP != NULL) {
			for (unsigned int i = 0; i < PConstP->size(); i++) {
				((*PConstP)[i])->Update(Order);
			}
		}

	}

	calcE();

	if (Order > 0) {
		calcdE();
	}
	if (Order > 1) {
		calcddE();
	}

}

void Rod::calcE() {

	calcEbend();
	calcEext();
	calcEgrav();

	E = Ebend + Egrav - Eext;

}

void Rod::calcEbend() {
	Ebend = 0.0;
	for (unsigned int i = 0; i < (*LCvecP).size(); i++) {
		Ebend += ((*LCvecP)[i]).getE();
	}
}

void Rod::calcEext() {
	Eext = 0.0;
	for (unsigned int i = 0; i < (*NvecP).size(); i++) {
		Eext += ((*NvecP)[i]).getX().dot(((*NvecP)[i]).getFnode());
	}
}

void Rod::calcEgrav() {
	Egrav = 0.0;
	for (unsigned int i = 0; i < (*NvecP).size(); i++) {
		Egrav += rho * ((*NvecP)[i]).getX().dot(-1.0 * gravdir);
	}

}

void Rod::calcdE() {
	if (Ndof == 0) {

		calcNdof();
		dE.resize(Ndof);
		dEbend.resize(Ndof);
		dEext.resize(Ndof);
		dEgrav.resize(Ndof);
	}

	calcdEbend();
	calcdEext();
	calcdEgrav();

	dE = dEbend + dEgrav - dEext;

}

void Rod::calcdEbend() {

	dEbend.setZero();

	for (unsigned int i = 0; i < (*LCvecP).size(); i++) {

		for (unsigned int j = 0; j < (((*LCvecP)[i]).getGlobalNodes()).size();
				j++) {

			dEbend((((*LCvecP)[i]).getGlobalNodes())[j]) +=
					(((*LCvecP)[i]).getdE())[j];
		}

	}
}

void Rod::calcdEext() {
	dEext.setZero();
	for (unsigned int i = 0; i < (*NvecP).size(); i++) {
		for (unsigned int j = 0; j < 3; j++) {
			dEext((((*NvecP)[i]).getGlobalNodes())[j]) +=
					(((*NvecP)[i]).getFnode())[j];
		}
	}

}

void Rod::calcdEgrav() {
	dEgrav.setZero();
	for (unsigned int i = 0; i < (*NvecP).size(); i++) {
		for (unsigned int j = 0; j < 3; j++) {
			dEgrav((((*NvecP)[i]).getGlobalNodes())[j]) -= (1.0) * rho
					* gravdir[j];
		}
	}

}

void Rod::calcddE() {

	//Eigen::SparseMatrix<double> ddE(Ndof,Ndof);

	calcddEbend();

	//ddE.setFromTriplets(ddEbendtrips.begin(),ddEbendtrips.end());
	ddEtrips = ddEbendtrips;

}

void Rod::calcddEbend() {
	//Eigen::SparseMatrix<double> ddEbend(Ndof,Ndof);
	ddEbend.setZero();

	//std::vector<T> ddEbendtrips;
	ddEbendtrips.clear();
	ddEbendtrips.reserve(121 * (*LCvecP).size());
	for (unsigned int i = 0; i < (*LCvecP).size(); i++) {
		int gnsize = (((*LCvecP)[i]).getGlobalNodes()).size();
		for (unsigned int j = 0; j < gnsize; j++) {
			for (unsigned int k = 0; k < gnsize; k++) {
				//int glj = (((*LCvecP)[i]).getGlobalNodes())[j];
				//int glk = (((*LCvecP)[i]).getGlobalNodes())[k];
				//ddEbend(glj,glk) += (((*LCvecP)[i]).getddE())(j,k);
				ddEbendtrips.push_back(
						T((((*LCvecP)[i]).getGlobalNodes())[j],
								(((*LCvecP)[i]).getGlobalNodes())[k],
								(((*LCvecP)[i]).getddE())(j, k)));
			}
		}
	}
	//ddEbend.setFromTriplets(ddEbendtrips.begin(),ddEbendtrips.end());

}
double Rod::getE() {
	return E;
}

double Rod::getEbend() {
	return Ebend;
}

double Rod::getEext() {
	return Eext;
}

double Rod::getEgrav() {
	return Egrav;
}

Eigen::VectorXd Rod::getdE() {
	return dE;
}

Eigen::VectorXd Rod::getdEbend() {
	return dEbend;
}

Eigen::VectorXd Rod::getdEext() {
	return dEext;
}

Eigen::VectorXd Rod::getdEgrav() {
	return dEgrav;
}

std::vector<T> Rod::getddE() {
	return ddEtrips;
}

Eigen::SparseMatrix<double> Rod::getddEbend() {
	return ddEbend;
}

void Rod::setrho(double r) {
	rho = r;
}

double Rod::getrho() {
	return rho;
}

void Rod::setNn(int nloc) {
	Nn = nloc;
}

int Rod::getNn() {
	return Nn;
}

void Rod::setNvecP(std::vector<Node> *np) {
	NvecP = np;
}

void Rod::setFvecP(std::vector<Frame> *fp) {
	FvecP = fp;
}

void Rod::setEvecP(std::vector<Element> *ep) {
	EvecP = ep;
}

void Rod::setMvecP(std::vector<Material> *mp) {
	MvecP = mp;
}

void Rod::setLCvecP(std::vector<LocalConstructor> *lcp) {
	LCvecP = lcp;
}

void Rod::setParallel(bool t) {
	Parallelized = t;
}
void Rod::setGravDir(Eigen::Vector3d &dir) {
	gravdir = dir.normalized();
}

void Rod::calcNdof() {
	std::vector<int> allGc;
	for (unsigned int i = 0; i < (*NvecP).size(); i++) {
		for (unsigned int j = 0; j < (((*NvecP)[i]).getGlobalNodes()).size();
				j++) {
			allGc.push_back((((*NvecP)[i]).getGlobalNodes())(j));

		}
	}
	/*
	 std::vector<int>::iterator ip;

	 int count;
	 sort(allGc.begin(), allGc.end());

	 // Using std::unique and std::distance to count
	 // unique elements in a container
	 count = std::distance(allGc.begin(),
	 std::unique(allGc.begin(), allGc.begin() + allGc.size()));
	 */
	int maxind = 0;
	for (unsigned int i = 0; i < allGc.size(); i++) {
		if (allGc[i] > maxind) {
			maxind = allGc[i];
		}
	}
	Ndof = maxind + 1;

}

void Rod::calcGlobaltoFree() {
	FreetoGlobal.clear();

	std::vector<pair<int, double*>> fnodes;
	for (unsigned int i = 0; i < (*NvecP).size(); i++) {
		if ((*NvecP)[i].isXFree()) {
			for (unsigned int j = 0; j < ((*NvecP)[i]).getGlobalNodes().size();
					j++) {
				FreetoGlobal.push_back((((*NvecP)[i]).getGlobalNodes())[j]);

			}

		}
	}
	for (unsigned int i = 0; i < (*FvecP).size(); i++) {
		if ((*FvecP)[i].isPhiFree()) {
			FreetoGlobal.push_back((((*FvecP)[i]).getGlobalNodes())[3]);

		}
	}

	std::sort(FreetoGlobal.begin(), FreetoGlobal.end());

	GlobaltoFree.clear();
	GlobaltoFree.resize(Ndof);

	for (int i = 0; i < GlobaltoFree.size(); i++) {
		std::pair<bool, int> find;
		find = findInVector(FreetoGlobal, i);
		GlobaltoFree[i] = find.second;
	}
	calcGenCoord();

}

std::vector<int> Rod::getGlobaltoFree() {
	return GlobaltoFree;
}
std::vector<int> Rod::getFreetoGlobal() {
	return FreetoGlobal;
}

template<typename T>
std::pair<bool, int> Rod::findInVector(const std::vector<T> &vecOfElements,
		const T &element) {
	std::pair<bool, int> result;
	// Find given element in vector
	auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);
	if (it != vecOfElements.end()) {
		result.second = distance(vecOfElements.begin(), it);
		result.first = true;
	} else {
		result.first = false;
		result.second = -1;
	}
	return result;
}
void Rod::calcGenCoord() {
	GenCoordAll.resize(Ndof);

	for (unsigned int i = 0; i < (*NvecP).size(); i++) {

		GenCoordAll[4 * i] = (*NvecP)[i].getX()[0];
		GenCoordAll[4 * i + 1] = (*NvecP)[i].getX()[1];
		GenCoordAll[4 * i + 2] = (*NvecP)[i].getX()[2];

		if (i < ((*NvecP).size() - 1)) {
			GenCoordAll[4 * i + 3] = (*FvecP)[i].getPhi();
		}
	}

	GenCoord.resize(FreetoGlobal.size());
	for (unsigned int i = 0; i < Ndof; i++) {
		if (GlobaltoFree[i] >= 0) {
			GenCoord[GlobaltoFree[i]] = GenCoordAll[i];
		}
	}

}
Eigen::VectorXd Rod::ExtractGenCoord() {
	calcGenCoord();
	return GenCoord;
}
Eigen::VectorXd Rod::ExtractGenCoordAll() {
	calcGenCoord();
	return GenCoordAll;
}
void Rod::ApplyGenCoordAll(Eigen::VectorXd &Xin) {
	for (unsigned int i = 0; i < (*NvecP).size(); i++) {
		Eigen::Vector3d xtemp;
		xtemp << Xin[4 * i], Xin[4 * i + 1], Xin[4 * i + 2];
		((*NvecP)[i]).setX(xtemp);
	}
	for (unsigned int i = 0; i < (*FvecP).size(); i++) {
		((*FvecP)[i]).setPhi(Xin[4 * i + 3]);
	}
}
void Rod::ApplyGenCoord(Eigen::VectorXd &Xin) {
	if (GenCoordAll.size() == 0) {
		calcGenCoord();
	}
	Eigen::VectorXd Xnew = GenCoordAll;
	for (unsigned int i = 0; i < Xin.size(); i++) {
		Xnew[FreetoGlobal[i]] = Xin[i];
	}
	ApplyGenCoordAll(Xnew);
}

void Rod::setdtoD() {
	for (unsigned int i = 0; i < (*FvecP).size(); i++) {
		((*FvecP)[i]).setD((*FvecP)[i].getd());
	}
}
void Rod::setPConstP(std::vector<Constraint*> *PC) {
	PConstP = PC;
	Nc = (*PConstP).size();
}
std::vector<Constraint*>* Rod::getPConstP() {
	return PConstP;
}
int Rod::getNc() {
	return Nc;
}
