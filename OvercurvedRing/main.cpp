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
#include <Constraint.h>
#include <Solver.h>
#include "read_file.h"
#include <fstream>
#include <iomanip>
#include "omp.h"
#include <gsl/gsl_sf_ellint.h>

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
void writevtk(std::string, std::vector<Node>&, std::vector<Frame>&);
void writelist(std::string, std::vector<double> &);

int main() {
	//std::vector<Node> nodevector(10);
	read_file rf;
	some_data dat;
	std::string fileNamestr("inputFile.in");

	char fileName[MAXLINE];
	//for (unsigned int i = 0; i < fileNamestr.length(); i++) {
	//	fileName[i] = fileNamestr[i];
	//}
	std::cout << "Please enter an input file: " << std::endl;
	std::cin >> fileName;
	//filename = { 'i', 'n', '' }
	rf.readInputFile(fileName, dat);

	int Nn = dat.Nn;
	int Nstep = dat.Nstep;
	//Eigen::Vector3d BMod = dat.BMod;
	double width = dat.width;
	double height = dat.height;
	double rho = dat.rho;
	Eigen::Vector3d fend = dat.fend;
	int SolverRead = dat.SolverTag;
	char *foldername = dat.FolderName;
	char *filename = dat.FileName;
	cout << foldername << endl;
	double pertsc = dat.pertsc;
	double OC = dat.phiend;
	

	Eigen::Vector3d BMod;
	BMod.setZero();
	//double Ix = pi*width*pow(height,3.0)/4.0;
	//double Iy = pi*height*pow(width,3.0)/4.0;
	//double J = Ix + Iy;
	
	double nu = 0.41;
	double Y = 0.256*height*pow(width,3.0);;
	double G = 1.0/(2.0*(1.0+nu));

	double Ix = width*pow(height,3.0)/12.0;
	double Iy = height*pow(width,3.0)/12.0;
	double J = G*Y;
	BMod << Ix, Iy, J;

	std::vector<Node> nvec(Nn);
	std::vector<double> Xlocs = Linspace(0.0,(double) Nn, Nn);
	std::vector<double> si = Linspace(0.0,1.0,nvec.size());
	for (unsigned int i = 0; i < nvec.size(); i++) {
		Eigen::Vector3d Xi;
		Xi << si[i], 0.0, 0.0;

		nvec[i].setX(Xi);
		Eigen::VectorXi gntemp(3);
		gntemp << (4 * i), (4 * i + 1), (4 * i + 2);
		nvec[i].setGlobalNodes(gntemp);
	}
	nvec[0].FreeNodeLoc(false);
	nvec[1].FreeNodeLoc(false);


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
		Eigen::VectorXi gntemp(7);
		gntemp << (4 * i), (4 * i + 1), (4 * i + 2), (4 * i + 3), (4 * i + 4), (4
			* i + 5), (4 * i + 6);
		fvec[i].setGlobalNodes(gntemp);
	}
	fvec[0].FreePhi(false);


	std::vector<Frame> *fvecptr;
	fvecptr = &fvec;


	std::vector<Element> evec(Nn - 2);
	for (unsigned int i = 0; i < evec.size(); i++) {
		evec[i].setNodesPointers(&nvec[i], &nvec[i + 1], &nvec[i + 2]);
		evec[i].setFramesPointers(&fvec[i], &fvec[i + 1]);
		VectorXi gctemp(11);
		gctemp << (4 * i), (4 * i + 1), (4 * i + 2), (4 * i + 3), (4 * i + 4), (4
			* i + 5), (4 * i + 6), (4 * i + 7), (4 * i + 8), (4 * i + 9), (4
			* i + 10);
			evec[i].setGlobalNodes(gctemp);
		}

		std::vector<Element> *evecptr;
		evecptr = &evec;

		std::vector<Material> mvec(evec.size());

		for (unsigned int i = 0; i < mvec.size(); i++) {
			Eigen::Vector3d BModtemp;
			BModtemp = BMod * pow(Nn - 2.0, 2.0);
			mvec[i].SetElementPointer(&evec[i]);
			mvec[i].setBMod(BModtemp);
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
		r1.setParallel(false);
		r1.Update(2);
		double rhosc = rho / ((double) Nn);
		r1.setrho(rhosc);
		r1.calcGlobaltoFree();

		std::vector<Inextensibility> Inext(fvecptr->size());

		for (unsigned int i = 0; i < Inext.size(); i++) {
			Inext[i].setFramePointer(&(fvec[i]));
			Inext[i].setdS(Xlocs[1] - Xlocs[0]);
			Inext[i].setGlobalNodes(fvec[i].getGlobalNodes());
		}
		std::vector<Constraint*> ConstP;
		std::vector<Constraint*> ConstP2;
		for (unsigned int i = 1; i < (Inext.size()); i++) {
			ConstP.push_back(&Inext[i]);
			if (i < (Inext.size() - 1)) {
				ConstP2.push_back(&Inext[i]);
			}
		}

		std::vector<Constraint*> *PConstP;
		std::vector<Constraint*> *PConstP2;
		PConstP = &ConstP;
		PConstP2 = &ConstP2;

		r1.setPConstP(PConstP);
		r1.Update(2);

		Solver s1;
		s1.setRodP(&r1);
		s1.setSolveDetailLevel(SolverRead);
	//s1.Solve();
		int Nstep1 = 50;
		std::vector<double> mult1 = Linspace(0.0,1.0,Nstep1);
		std::vector<double> mult = Linspace(0.0, 1.0, Nstep);
		cout << "Bending Modulus: " << BMod[0] << ", " << BMod[1] << ", " << BMod[2] << endl; 
//Closing Circle
	for (unsigned int i = 0; i < Nstep1; i++) {
		cout << "Step: " << i << endl;

		Eigen::Vector3d k0temp;
		k0temp.setZero();
		k0temp[0] = mult1[i]*2.0*pi/Nn;
			for (unsigned int j = 0; j < mvec.size(); j++){
			mvec[j].setkappa0(k0temp);
		}
		//cout << mult[i] << endl;
		double rhotemp = mult1[i] * rhosc;
		//nvec.back().setFnode(mult[i] * fend);


		r1.setrho(0.0*rhotemp);
		r1.Update(2);
		s1.Solve();
		r1.setdtoD();
		//cout << "Energy: " << r1.getE() / Nn << endl;
		//cout << "Last Node Pos: " << nvec.back().getX()[0] / (Nn - 1) << " "
		//		<< nvec.back().getX()[1] / (Nn - 1) << " "
		//		<< nvec.back().getX()[2] / (Nn - 1) << endl;

		//cout << r1.getdEbend() << endl;

		std::string directoryname;
		std::ostringstream ss;
		ss << foldername << "/" << filename << std::setw(5) << std::setfill('0')
		<< i << ".vtk";

		writevtk(ss.str(), nvec, fvec);
	}


	fvec.back().FreePhi(false);
	nvec.back().FreeNodeLoc(false);
	nvec[nvec.size()-2].FreeNodeLoc(false);
	Eigen::Vector3d X0 = nvec[0].getX();
	Eigen::Vector3d X0p = nvec[1].getX();

	nvec.back().setX(X0p);
	nvec[nvec.size()-2].setX(X0);


	r1.setPConstP(PConstP2);
	r1.calcGlobaltoFree();
	Solver s2;
	s2.setRodP(&r1);
	s2.setSolveDetailLevel(SolverRead);
	Eigen::Vector3d k0temp;
	double k0max1 = OC*2.0*pi/((double) Nn);
	std::vector<double> k0vals1 = Linspace(2.0*pi/((double) Nn), k0max1, Nstep); 

	//Applying gravity to introduce asymmetry
	k0temp.setZero();
		k0temp[0] = 2.0*pi/((double) Nn);
			
		for (unsigned int j = 0; j < mvec.size(); j++){
			mvec[j].setkappa0(k0temp);
		}

	for (unsigned int i = 0; i < Nstep; i++) {
		cout << "Step: " << i << endl;
			

		//cout << mult[i] << endl;
		double rhotemp = mult[i] * rhosc;
		//nvec.back().setFnode(mult[i] * fend);


		r1.setrho(rhotemp);
		r1.Update(2);
		s2.Solve();
		r1.setdtoD();

		std::string directoryname;
		std::ostringstream ss;
		ss << foldername << "/" << filename << std::setw(5) << std::setfill('0')
			<< i << ".vtk";

		writevtk(ss.str(), nvec, fvec);
	}
	// Also increasing curvature

	for (unsigned int i = 0; i < Nstep; i++) {
		cout << "Step: " << i << endl;
			
		k0temp.setZero();
		cout << k0vals1[i] << endl;
		k0temp[0] = k0vals1[i];
			
		for (unsigned int j = 0; j < mvec.size(); j++){
			mvec[j].setkappa0(k0temp);
		}

		//nvec.back().setFnode(mult[i] * fend);



		r1.Update(2);
		s2.Solve();
		r1.setdtoD();

		std::string directoryname;
		std::ostringstream ss;
		ss << foldername << "/" << filename << std::setw(5) << std::setfill('0')
			<< i << ".vtk";

		writevtk(ss.str(), nvec, fvec);
	}

// removing gravity
	for (unsigned int i = 0; i < Nstep; i++) {
		cout << "Step: " << Nstep + i << endl;


		//cout << mult[i] << endl;
		double rhotemp = (1.0-mult[i]) * rhosc;
		//nvec.back().setFnode(mult[i] * fend);


		r1.setrho(rhotemp);
		r1.Update(2);
		s2.Solve();
		r1.setdtoD();

		std::string directoryname;
		std::ostringstream ss;
		ss << foldername << "/" << filename << std::setw(5) << std::setfill('0')
			<< Nstep + i << ".vtk";

		writevtk(ss.str(), nvec, fvec);
	}



	// saving the current twist to be add as an offset
	std::vector<double> k2vals(evec.size());
	std::vector<double> k3vals(evec.size());
	for (unsigned int i = 0; i < evec.size(); i++) {
		Eigen::Vector3d kvalstemp = evec[i].getkappa();
		k2vals[i] = kvalstemp[1];
		k3vals[i] = kvalstemp[2];
	}
	double k0sc = .01;
		



//lowering the overcurve
		std::vector<double> Ovals(Nstep);
		std::vector<double> Dvals(Nstep);
		std::vector<double> Hvals(Nstep);




		double k0max = k0temp[0];
		double k0O1 = 2.0*pi/((double) Nn);
		std::vector<double> k0vals = Linspace(k0max,k0O1,Nstep);


		for (unsigned int i = 0; i < Nstep; i++){
			Ovals[i] = k0vals[i]/k0O1;
		}
				double Dmin = 2.0*((double) Nn);
		int indD1 = 0;
		int Doffset = nearbyint(Nn/2);
		int indD2 = indD1 + Doffset;
		Eigen::Vector3d Dnvec;

	for (unsigned int i = 0; i < Nstep; i++) {
		cout << "Step: " << 2*Nstep - i-1 << endl;


		//cout << mult[i] << endl;

		//nvec.back().setFnode(mult[i] * fend);
		k0temp.setZero();
		k0temp[0] = k0vals[i];
			
		for (unsigned int j = 0; j < mvec.size(); j++){
			k0temp[1] = k0sc*k2vals[j];
			k0temp[2] = k0sc*k3vals[j];
			mvec[j].setkappa0(k0temp);
		}

		r1.setrho(0.0);
		r1.Update(2);
		s2.Solve();
		r1.setdtoD();



		//Calculating the D value


		if (i == 0){
			for (unsigned int j = 0; j < Doffset; j++){
				Eigen::Vector3d Dvec = nvec[j+Doffset].getX() - nvec[j].getX();
				double Dtemp = Dvec.norm();
				Dnvec = Dvec/Dtemp;
				if (Dtemp < Dmin) {
					indD1 = j;
					indD2 = j+ Doffset;
					Dmin = Dtemp;
				}
			}
		} else {
			Eigen::Vector3d Dvec = nvec[indD2].getX() - nvec[indD1].getX();
			Dmin = Dvec.norm();

		}
		cout << indD1 << endl;
		Dvals[i] = Dmin;


		std::string directoryname;
		std::ostringstream ss;
		ss << foldername << "/" << filename << std::setw(5) << std::setfill('0')
			<< 2*Nstep - i -1 << ".vtk";

		writevtk(ss.str(), nvec, fvec);

	}
	writelist("Dlist.txt",Dvals);
	writelist("Olist.txt",Ovals);



	for (unsigned int i = 0; i < Nn; i++) {
		cout << nvec[i].getX()[0] / (Nn - 1) << " "
			<< nvec[i].getX()[1] / (Nn - 1) << " "
			<< nvec[i].getX()[2] / (Nn - 1) << endl;
	}

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
		std::vector<Frame> &fvec) {
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

		myfile.close();

	}

	void writelist(std::string filename, std::vector<double> &list) {
		ofstream myfile;
		myfile.open(filename);
		

	// Writing dataset info list

		for (unsigned int i = 0; i < list.size(); i++) {
			myfile << list[i]  << "\n";
		}

		myfile.close();

	}

