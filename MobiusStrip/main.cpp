#include <iostream>
#include <quaternion.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <math.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Frame.h>
#include <Element.h>
#include <Material.h>
#include <LocalConstructor.h>
#include <Rod.h>
#include <Inextensibility.h>
#include <k20.h>
#include <RibbonConst.h>
#include <Constraint.h>
#include <Solver.h>
#include <EtapLimP.h>
#include <EtapLimN.h>
#include "read_file.h"
#include <fstream>
#include <iomanip>
#include "omp.h"
#include <X3Const.h>

using namespace std;
using namespace Eigen;
typedef Eigen::Triplet<double> T;
#include <chrono>
using namespace std::chrono;
double pi = 3.1415926535897932384626433;
void Print3D(const Eigen::Tensor<double, 5>&, int, int);
void Print3D2(const Eigen::Tensor<double, 3>&);
void PrintTriplet(T);
Eigen::MatrixXd TripList2EMat(std::vector<T>);
std::vector<double> Linspace(double, double, int);
void writevtk(std::string, std::vector<Node>&, std::vector<Frame>&, double);

int main() {
	//std::vector<Node> nodevector(10);
	read_file rf;
	some_data dat;
	std::string fileNamestr("inputFile.in");
	//char fileName[MAXLINE];
	//for (unsigned int i = 0; i < fileNamestr.length(); i++) {
	//	fileName[i] = fileNamestr[i];
	//}
	//std::cout << "Please enter an input file: " << std::endl;
	//std::cin >> fileName;
	char fileName[] = "inputFile.in";
	//filename = { 'i', 'n', '' }
	rf.readInputFile(fileName, dat);

	int Nn = dat.Nn;
	int Nstep = dat.Nstep;
	Eigen::VectorXd BMod = dat.BMod;
	double rho = dat.rho / ((double) Nn);
	Eigen::Vector3d fend = dat.fend;
	int SolverRead = dat.SolverTag;
	char *foldername = dat.FolderName;
	char *filename = dat.FileName;
	cout << foldername << endl;
	double width = dat.width * Nn;
	int NumThreads = dat.NumThreads;

	double T = dat.T;
	std::vector<double> tsteps = Linspace(0.0, T, Nstep);
	double dt = tsteps[1] - tsteps[0];
	double mu = dat.viscosity;

	double Light_Intensity = dat.LightIntensity / ((double) Nn);
	Eigen::Vector3d IllumDirec;
	IllumDirec = dat.LightDirec;
	IllumDirec.normalize();

	double w10 = dat.w0;
	double eta0 = dat.eta0;
	Eigen::Matrix3d Kr0;
	Kr0.setZero();

	Eigen::Matrix3d Dk;
	Dk << -w10 * pow(eta0, 2.0), 0.0, w10 * eta0, 0.0, 0.0, 0.0, w10 * eta0, 0.0, -w10;
	Dk = Dk / sqrt(pow(eta0, 2) + 1.0);
	Kr0 = -Dk / ((double) Nn - 1.0);
	std::vector<Node> nvec(Nn);
	std::vector<double> Xvals = Linspace(0.0, 1.0 * Nn, Nn);
	for (unsigned int i = 0; i < nvec.size(); i++) {
		Eigen::Vector3d Xi;
		Xi << Xvals[i], 0.0, 0.0;

		nvec[i].setX(Xi);
		Eigen::VectorXi gntemp(4);
		gntemp << (5 * i), (5 * i + 1), (5 * i + 2), (5 * i + 3);
		nvec[i].setGlobalNodes(gntemp);
	}
	nvec[0].FreeNodeLoc(false);
	nvec[1].FreeNodeLoc(false);
	nvec[0].FreeEtaSet(false);
	nvec[1].FreeEtaSet(false);
	//nvec[1].FreeEtaSet(false);
	nvec.back().FreeEtaSet(false);

	std::vector<Node> *nvecptr;
	nvecptr = &nvec;

	std::vector<Frame> fvec(Nn - 1);

	Eigen::VectorXd phivals(Nn - 1);
	phivals.setZero();

	for (unsigned int i = 0; i < fvec.size(); i++) {
		fvec[i].setNodesPointers(&nvec[i], &nvec[i + 1]);
		Quat<double> qtemp;
		qtemp.Set(0.5, 0.5, 0.5, 0.5);
		fvec[i].setD(qtemp);
		fvec[i].setPhi(phivals[i]);
		Eigen::VectorXi gntemp(9);
		gntemp << (5 * i), (5 * i + 1), (5 * i + 2), (5 * i + 3), (5 * i + 4), (5
				* i + 5), (5 * i + 6), (5 * i + 7), (5 * i + 8);
		fvec[i].setGlobalNodes(gntemp);
	}
	fvec[0].FreePhi(false);

	std::vector<Frame> *fvecptr;
	fvecptr = &fvec;

	std::vector<Element> evec(Nn - 2);
	for (unsigned int i = 0; i < evec.size(); i++) {
		evec[i].setNodesPointers(&nvec[i], &nvec[i + 1], &nvec[i + 2]);
		evec[i].setFramesPointers(&fvec[i], &fvec[i + 1]);
		VectorXi gctemp(14);
		gctemp << (5 * i), (5 * i + 1), (5 * i + 2), (5 * i + 3), (5 * i + 4), (5
				* i + 5), (5 * i + 6), (5 * i + 7), (5 * i + 8), (5 * i + 9), (5
				* i + 10), (5 * i + 11), (5 * i + 12), (5 * i + 13);
		evec[i].setGlobalNodes(gctemp);
	}

	std::vector<Element> *evecptr;
	evecptr = &evec;

	std::vector<Material> mvec(evec.size());

	for (unsigned int i = 0; i < mvec.size(); i++) {
		Eigen::VectorXd BModtemp;
		BModtemp = BMod * pow(Nn - 2.0, 2.0) / width;
		mvec[i].SetElementPointer(&evec[i]);
		mvec[i].setBMod(BModtemp);
		mvec[i].setw(width);
		mvec[i].setnu(0.0);

	}

	std::vector<Material> *mvecptr;
	mvecptr = &mvec;

	std::vector<LocalConstructor> lcvec(evec.size());

	for (unsigned int i = 0; i < lcvec.size(); i++) {
		lcvec[i].setElementPointer(&evec[i]);
		lcvec[i].setMaterialPointer(&mvec[i]);
		lcvec[i].setGlobalNodes(evec[i].getGlobalNodes());
	}

	std::vector<LocalConstructor> *lcvecptr;
	lcvecptr = &lcvec;

	Rod r1;
	r1.setNvecP(nvecptr);
	r1.setFvecP(fvecptr);
	r1.setEvecP(evecptr);
	r1.setMvecP(mvecptr);
	r1.setLCvecP(lcvecptr);
	if (NumThreads == 0) {
		r1.setParallel(false);
	} else {
		r1.setParallel(true);
		r1.SetNumThreads(NumThreads);
	}

	r1.Update(2);
	double rhosc = rho / ((double) Nn);
	r1.setrho(rhosc);
	r1.calcGlobaltoFree();

	std::vector<Inextensibility> Inext(fvecptr->size());

	std::vector<X3Const> XendConst(nvecptr->size());
	for (unsigned int i = 0; i < XendConst.size(); i++) {
		XendConst[i].setNodePointer(&nvec[i]);
		XendConst[i].setX3min(-1.0 * ((double) Nn) / 3.0);
		XendConst[i].setGlobalNodes(nvec[i].getGlobalNodes());
	}

	std::vector<EtapLimP> EtapLimPConst(mvec.size());
	std::vector<EtapLimN> EtapLimNConst(mvec.size());
	double BarrierLoc = 1.99;
	for (unsigned int i = 0; i < mvec.size(); i++){
		EtapLimPConst[i].setMaterialPointer(&mvec[i]);
		EtapLimNConst[i].setMaterialPointer(&mvec[i]);

		EtapLimPConst[i].setBarrierLoc(BarrierLoc);
		EtapLimNConst[i].setBarrierLoc(BarrierLoc);

		Eigen::VectorXi GNloc(3);
		GNloc[0] = (evec[i].getGlobalNodes())[3];
		GNloc[1] = (evec[i].getGlobalNodes())[8];
		GNloc[2] = (evec[i].getGlobalNodes())[13];

		EtapLimPConst[i].setGlobalNodes(GNloc);
		EtapLimNConst[i].setGlobalNodes(GNloc);

	}


	for (unsigned int i = 0; i < Inext.size(); i++) {
		Inext[i].setFramePointer(&(fvec[i]));
		Inext[i].setdS(1.0);
		Inext[i].setGlobalNodes(fvec[i].getGlobalNodes());
	}

	std::vector<k20> k20consts(evecptr->size());
	for (unsigned int i = 0; i < k20consts.size(); i++) {
		k20consts[i].setElementPointer(&(evec[i]));
		k20consts[i].setdS(1.0);
		k20consts[i].setGlobalNodes(evec[i].getGlobalNodes());
	}
	std::vector<RibbonConst> RibbonConsts(evecptr->size());
	for (unsigned int i = 0; i < k20consts.size(); i++) {
		RibbonConsts[i].setElementPointer(&(evec[i]));
		RibbonConsts[i].setdS(1.0);
		RibbonConsts[i].setGlobalNodes(evec[i].getGlobalNodes());
	}
	std::vector<Constraint*> ConstP;
	std::vector<Constraint*> ConstP2;
	std::vector<Constraint*> ConstInP;

	for (unsigned int i = 1; i < Inext.size(); i++) {
		ConstP.push_back(&Inext[i]);
		if (i < (Inext.size() - 1)) {
			ConstP2.push_back(&Inext[i]);
		}
	}
	for (unsigned int i = 0; i < k20consts.size(); i++) {
		ConstP.push_back(&k20consts[i]);
		ConstP2.push_back(&k20consts[i]);

	}
	for (unsigned int i = 0; i < RibbonConsts.size(); i++) {
		ConstP.push_back(&RibbonConsts[i]);
		ConstP2.push_back(&RibbonConsts[i]);
	}

	for (unsigned int i = 0; i < EtapLimPConst.size(); i++){
		ConstInP.push_back(&EtapLimPConst[i]);
		ConstInP.push_back(&EtapLimNConst[i]);
	}

	std::vector<Constraint*> *PConstP;
	std::vector<Constraint*> *PConstP2;
	std::vector<Constraint*> *PConstInP;
	PConstP = &ConstP;
	PConstP2 = &ConstP2;
	PConstInP = &ConstInP;

	r1.setPConstP(PConstP);
	r1.setPConstInP(PConstInP);
	r1.Update(2);

	Solver s1;
	s1.setRodP(&r1);
	s1.setSolveDetailLevel(SolverRead);
	//s1.Solve();
	int Nstep1 = 50;
	int Nstep2 = 500;

	std::vector<double> mult = Linspace(0.0, 1.0, Nstep);
	std::vector<double> mult1 = Linspace(0.0, 1.0, Nstep1);
	std::vector<double> mult2 = Linspace(0.0,1.0, Nstep2);

	r1.SetXtoXold();
	double dphi = .5 * pi / ((double) Nstep2);

	cout << "dphi: " << dphi << endl;
	r1.calcGlobaltoFree();
	for (unsigned int i = 0; i < Nstep1; i++) {
		cout << "Step: " << i << endl;
		//double rhotemp = mult1[i] * rhosc * width;
		double rhotemp = rhosc * width;
		Eigen::Matrix3d K0temp = mult1[i] * Kr0;
		for (unsigned int j = 0; j < mvec.size(); j++) {
			mvec[j].setKr(K0temp);
		}

		fvec.size();

		r1.setrho(rhotemp);

		r1.Update(2);
		s1.Solve();
		r1.setdtoD();
		r1.SetXtoXold();
		//for (unsigned int j = 0; j < fvec.size(); j++) {
		//	cout << "Phi " << j << ": " << fvec[j].isPhiFree() << ", "
		//			<< fvec[j].getPhi() << endl;
		//}

		std::string directoryname;
		std::ostringstream ss;
		ss << foldername << "/" << filename << std::setw(5) << std::setfill('0')
				<< i << ".vtk";
		writevtk(ss.str(), nvec, fvec, dat.width);
	}

	// want to twist the end
	fvec.back().FreePhi(false);
	nvec.back().FreeNodeLoc(false);
	int Nendm1 = Nn - 2;
	nvec[Nendm1].FreeNodeLoc(false);
	nvec.back().FreeEtaSet(false);
	Eigen::Vector3d X0 = nvec[0].getX();
	nvec.back().setX(X0);
	Eigen::Vector3d dX = nvec[1].getX() - nvec[0].getX();
	Eigen::Vector3d Xnew = nvec[0].getX() - dX;
	nvec[Nendm1].setX(Xnew);
	r1.calcGlobaltoFree();
	r1.setPConstP(PConstP2);

	//double dphi = 2.0 * pi / ((double) Nstep);
	r1.SetNodeViscosity(mu * dt / (pow((double) Nn, 2.0)));
	Solver s2;
	s2.setRodP(&r1);
	s2.setSolveDetailLevel(SolverRead);
	Eigen::Matrix3d K0temp = 0.0 * Kr0;
	for (unsigned int i = 0; i < mvec.size(); i++) {
		mvec[i].setKr(K0temp);
	}

	for (unsigned int i = 0; i < Nstep2; i++) {
		cout << "Step: " << Nstep2 + i << endl;

		fvec.back().setPhi(dphi);
		fvec[0].setPhi(-1.0 * dphi);
		r1.Update(2);
		s2.Solve();
		r1.setdtoD();
		r1.SetXtoXold();

		std::string directoryname;
		std::ostringstream ss;
		ss << foldername << "/" << filename << std::setw(5) << std::setfill('0')
				<< Nstep2 + i << ".vtk";
		writevtk(ss.str(), nvec, fvec, dat.width);
	}
	fvec.back().setPhi(0.0);
	fvec[0].setPhi(0.0);
	r1.SetNodeViscosity(0.0 * mu * dt / (pow((double) Nn, 2.0)));

	/*
	 for (unsigned int i = 0; i < Nstep; i++) {
	 cout << "Step: " << 2 * Nstep + i << endl;

	 r1.Update(2);
	 s2.Solve();
	 r1.setdtoD();
	 r1.SetXtoXold();

	 std::string directoryname;
	 std::ostringstream ss;
	 ss << foldername << "/" << filename << std::setw(5) << std::setfill('0')
	 << 2 * Nstep + i << ".vtk";
	 writevtk(ss.str(), nvec, fvec, dat.width);
	 }
	 */

	double wmax = 1.0 / pi;
	fvec.back().setPhi(0.0);
	fvec[0].setPhi(0.0);
	std::vector<double> wvals = Linspace(dat.width, wmax, Nstep);

		for (unsigned int j = 0; j < mvec.size(); j++) {
			//mvec[j].setw(((double) Nn) * wvals[i]);
			Eigen::VectorXd BModtemp;
			BModtemp = BMod * pow(Nn - 2.0, 2.0) / width;;
			BModtemp[3] = 0.0;
			mvec[j].setBMod(BModtemp);
		}
	for (unsigned int i = 0; i < Nstep; i++) {
		cout << "Step: " << 2 * Nstep + i << endl;



		for (unsigned int j = 0; j < mvec.size(); j++) {
			mvec[j].setw(((double) Nn) * wvals[i]);
			//Eigen::VectorXd BModtemp;
			//BModtemp = BMod * pow(Nn - 2.0, 2.0)/(((double) Nn) * wvals[i]);
			//BModtemp[3] *= wvals[i];
			//mvec[j].setBMod(BModtemp);
		}


		r1.Update(2);
		s2.Solve();
		r1.setdtoD();
		r1.SetXtoXold();

		std::string directoryname;
		std::ostringstream ss;
		ss << foldername << "/" << filename << std::setw(5) << std::setfill('0')
				<< 2 * Nstep + i << ".vtk";
		writevtk(ss.str(), nvec, fvec, wvals[i]);

		ofstream k1;
		std::ostringstream k1n;
		k1n << "w_k1/" << "w_" << wvals[i] << "_k1.txt";
		k1.open(k1n.str());

		for (unsigned int i = 0; i < evec.size(); i++) {
			k1 << (float) evec[i].getkappa()[0] << "\n";
		}

		k1.close();

		ofstream k3;
		std::ostringstream k3n;
		k3n << "w_k3/" << "w_" << wvals[i] << "_k3.txt";
		k3.open(k3n.str());
		for (unsigned int i = 0; i < evec.size(); i++) {
			k3 << (float) evec[i].getkappa()[2] << "\n";
		}

		k3.close();

		ofstream eta;
		std::ostringstream etan;
		etan << "w_eta/" << "w_" << wvals[i] << "_eta.txt";
		eta.open(etan.str());
		for (unsigned int i = 0; i < nvec.size(); i++) {
			eta << (float) nvec[i].geteta() << "\n";
		}
		eta.close();

	}

	for (unsigned int i = 0; i < Nn; i++) {
		cout << nvec[i].getX()[0] / (Nn - 1) << " "
				<< nvec[i].getX()[1] / (Nn - 1) << " "
				<< nvec[i].getX()[2] / (Nn - 1) << endl;
	}

	// Saving curvature data after code runs

	ofstream k1;
	k1.open("k1.txt");

	for (unsigned int i = 0; i < evec.size(); i++) {
		k1 << (float) evec[i].getkappa()[0] << "\n";
	}

	k1.close();

	ofstream k3;
	k3.open("k3.txt");
	for (unsigned int i = 0; i < evec.size(); i++) {
		k3 << (float) evec[i].getkappa()[2] << "\n";
	}

	k3.close();

	ofstream eta;
	eta.open("eta.txt");
	for (unsigned int i = 0; i < nvec.size(); i++) {
		eta << (float) nvec[i].geteta() << "\n";
	}
	eta.close();

	return 0;
}

Eigen::MatrixXd TripList2EMat(std::vector<T> tlist) {
	int N = 0;
	int M = 0;
	for (unsigned int i = 0; i < tlist.size(); i++) {
		if (tlist[i].row() > N) {
			N = tlist[i].row();

		}
		if (tlist[i].col() > M) {
			M = tlist[i].col();
		}
	}
	Eigen::MatrixXd Matout(N + 1, M + 1);
	Matout.setZero();
	for (unsigned int i = 0; i < tlist.size(); i++) {
		Matout(tlist[i].row(), tlist[i].col()) += tlist[i].value();
	}
	return Matout;
}
void PrintTriplet(T trip) {
	cout << trip.row() << ", " << trip.col() << ", " << trip.value() << endl;
}
void Print3D(const Eigen::Tensor<double, 5> &M, int a, int b) {
	for (unsigned int k = 0; k < 3; k++) {
		for (unsigned int j = 0; j < 3; j++) {
			cout << M(j, 0, k, a, b) << " " << M(j, 1, k, a, b) << " "
					<< M(j, 2, k, a, b) << endl;
		}
		cout << "---------------" << endl;
	}
}

void Print3D2(const Eigen::Tensor<double, 3> &M) {
	for (unsigned int k = 0; k < 2; k++) {
		for (unsigned int j = 0; j < 3; j++) {
			cout << M(j, 0, k) << " " << M(j, 1, k) << endl;
		}
		cout << "---------------" << endl;
	}
}

std::vector<double> Linspace(double a, double b, int N) {
	std::vector<double> vout(N);
	double Ndoub = 1.0 * N;
	for (unsigned int i = 0; i < N; i++) {
		vout[i] = a * ((Ndoub - 1.0) - i) / (Ndoub - 1.0)
				+ b * i / (Ndoub - 1.0);
	}
	return vout;
}

void writevtk(std::string filename, std::vector<Node> &nvec,
		std::vector<Frame> &fvec, double w) {
	ofstream myfile;
	myfile.open(filename);
	myfile << "# vtk DataFile Version 2.0\n";
	myfile << "Written badly by Kevin Korner\n";
	myfile << "ASCII\n";

	// Writing dataset info (X,Y,Z)

	myfile << "DATASET POLYDATA\n";
	myfile << "POINTS " << nvec.size() << " FLOAT\n";
	for (unsigned int i = 0; i < nvec.size(); i++) {
		myfile << nvec[i].getX()[0] / (nvec.size()) << " "
				<< nvec[i].getX()[1] / (nvec.size()) << " "
				<< nvec[i].getX()[2] / (nvec.size()) << "\n";
	}

	myfile << "\n";
	// Writing the connectivity
	myfile << "Lines " << nvec.size() - 1 << " " << 3 * (nvec.size() - 1)
			<< "\n";

	for (unsigned int i = 0; i < (nvec.size() - 1); i++) {
		myfile << 2 << " " << i << " " << i + 1 << "\n";
	}

	// Writing the director fields
	//d1
	myfile << "\n";

	myfile << "POINT_DATA " << nvec.size() << "\n";
	myfile << "VECTORS" << " d1 " << "FLOAT\n";
	Quat<double> e1quat;
	e1quat.Set(0.0, 1.0, 0.0, 0.0);
	Eigen::Vector3d d1temp;
	for (unsigned int i = 0; i < (fvec.size()); i++) {
		d1temp = (fvec[i].getd() * e1quat * fvec[i].getd().conj()).Getv();
		myfile << d1temp[0] << " " << d1temp[1] << " " << d1temp[2] << "\n";
	}
	myfile << d1temp[0] << " " << d1temp[1] << " " << d1temp[2] << "\n";

	//d2
	myfile << "\n";

	myfile << "VECTORS" << " d2 " << "FLOAT\n";
	Quat<double> e2quat;
	e2quat.Set(0.0, 0.0, 1.0, 0.0);
	Eigen::Vector3d d2temp;
	for (unsigned int i = 0; i < (fvec.size()); i++) {
		d2temp = (fvec[i].getd() * e2quat * fvec[i].getd().conj()).Getv();
		myfile << d2temp[0] << " " << d2temp[1] << " " << d2temp[2] << "\n";
	}
	myfile << d2temp[0] << " " << d2temp[1] << " " << d2temp[2] << "\n";

	//generatrix
	myfile << "\n";

	myfile << "VECTORS" << " eta " << "FLOAT\n";
	Quat<double> e3quat;
	e3quat.Set(0.0, 0.0, 0.0, 1.0);
	Eigen::Vector3d d3temp;
	Eigen::Vector3d etatemp;
	for (unsigned int i = 0; i < (fvec.size()); i++) {
		d1temp = (fvec[i].getd() * e1quat * fvec[i].getd().conj()).Getv();
		d3temp = (fvec[i].getd() * e3quat * fvec[i].getd().conj()).Getv();
		etatemp = (w) * (d1temp + nvec[i].geteta() * d3temp);
		myfile << etatemp[0] << " " << etatemp[1] << " " << etatemp[2] << "\n";
	}
	etatemp = (w) * (d1temp + nvec.back().geteta() * d3temp);
	myfile << etatemp[0] << " " << etatemp[1] << " " << etatemp[2] << "\n";

	myfile.close();

}

