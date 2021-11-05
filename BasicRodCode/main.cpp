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
	// In this code, we set up a basic running example of the code. It is based off the simulations
	// of the overcurvature, but will only include a basic simulation with curvature and end forces

	// Reading data in from the input file. This input file is there to allow for
	// easily changing certain parameters without having to recompile the code
	read_file rf;
	some_data dat;


	char fileName[MAXLINE];
	// Reads in the input file name from the command line. In our case it is inputFile.in
	std::cout << "Please enter an input file: " << std::endl;
	std::cin >> fileName;

	rf.readInputFile(fileName, dat);


	// Saving data from the input file to local variables for use
	int Nn = dat.Nn; // Number of nodes
	int Nstep = dat.Nstep; // Number of steps in iterative solves
	double width = dat.width; // Width used for bending moduli
	double height = dat.height; // Height used for bending moduli
	double L = dat.length;
	double rho = dat.rho; // Gravitationally applied load
	double YoungsMod = dat.youngsmod;
	Eigen::Vector3d fend = dat.fend; // Force applied to end
	int SolverRead = dat.SolverTag;
	char *foldername = dat.FolderName; // Folder name to save data to
	char *filename = dat.FileName; // File name for save data
	cout << foldername << endl;

	// Setting up a bending modulus vector. The implementation is that of a Kirchhoff Beam
	// the values are those given in the overcurvature example
	Eigen::Vector3d BMod; // initializing a constant Bending Modulus vector
	BMod.setZero();


	double nu = 0.41;
	double Y = 0.256*height*pow(width,3.0);
	double G = 1.0/(2.0*(1.0+nu));


	double Ix = YoungsMod*width*pow(height,3.0)/12.0;
	double Iy = YoungsMod*height*pow(width,3.0)/12.0;
	double J = YoungsMod*G*Y;
	BMod << Ix, Iy, J;

	// This section of code sets up many of the data structures necessary for code function.
	// Nodes, frames, elements, etc. are stored as std::vectors and their relations are stored
	std::vector<Node> nvec(Nn); // Vector of node objects
	std::vector<double> si = Linspace(0.0,L,nvec.size()); //Arc length parameters for each node
	for (unsigned int i = 0; i < nvec.size(); i++) {
		Eigen::Vector3d Xi;
		Xi << si[i], 0.0, 0.0; // setting initial positions. The rod is aligned with the x axis

		nvec[i].setX(Xi);
		Eigen::VectorXi gntemp(3);
		gntemp << (4 * i), (4 * i + 1), (4 * i + 2); // Setting the global degrees of freedom
		nvec[i].setGlobalNodes(gntemp);
	}
	nvec[0].FreeNodeLoc(false); //Fixes the first node in place (x1,x2,x3) are constant
	nvec[1].FreeNodeLoc(false); //Fixes the second node in place (x1,x2,x3) are constant


	std::vector<Node> *nvecptr;
	nvecptr = &nvec; // Node vector manipulated using a pointer to the std::vector

	std::vector<Frame> fvec(Nn - 1); //Initializing the frame vector

	Eigen::VectorXd phivals(Nn - 1); //Initializing the initial phi values
	phivals.setZero();

	for (unsigned int i = 0; i < fvec.size(); i++) {
		fvec[i].setNodesPointers(&nvec[i], &nvec[i + 1]); //linking the frame to the nodes
		Quat<double> qtemp;
		qtemp.Set(0.5, 0.5, 0.5, 0.5);
		fvec[i].setD(qtemp); // setting the initial D quaternion
		fvec[i].setPhi(phivals[i]); // setting the initial values of phi
		Eigen::VectorXi gntemp(7);
		gntemp << (4 * i), (4 * i + 1), (4 * i + 2), (4 * i + 3), (4 * i + 4), (4
				* i + 5), (4 * i + 6);
		fvec[i].setGlobalNodes(gntemp); // Setting the global degress of freedom associated with the frame
	}
	fvec[0].FreePhi(false); // Fixes the first frame to not rotate


	std::vector<Frame> *fvecptr;
	fvecptr = &fvec;


	std::vector<Element> evec(Nn - 2); // initializing the element, comprising of two frames
	for (unsigned int i = 0; i < evec.size(); i++) {
		evec[i].setNodesPointers(&nvec[i], &nvec[i + 1], &nvec[i + 2]); // setting the associated nodes
 		evec[i].setFramesPointers(&fvec[i], &fvec[i + 1]); // setting the associated frames
		VectorXi gctemp(11);
		gctemp << (4 * i), (4 * i + 1), (4 * i + 2), (4 * i + 3), (4 * i + 4), (4
				* i + 5), (4 * i + 6), (4 * i + 7), (4 * i + 8), (4 * i + 9), (4
						* i + 10); //setting the global degrees of freedom
		evec[i].setGlobalNodes(gctemp);
	}

	std::vector<Element> *evecptr;
	evecptr = &evec;

	std::vector<Material> mvec(evec.size()); //initializing the material vector.

	for (unsigned int i = 0; i < mvec.size(); i++) {
		Eigen::Vector3d BModtemp;
		BModtemp = BMod;
		mvec[i].SetElementPointer(&evec[i]); // Associating the material element with the element
		mvec[i].setBMod(BModtemp); // Setting the bending modulus. in this case constant
		double ds = si[i+1] - si[i];
		mvec[i].setli(ds); // Feeding in the node separation
	}

	std::vector<Material> *mvecptr;
	mvecptr = &mvec;

	std::vector<LocalConstructor> lcvec(evec.size()); // Local constructor calculates local energies and gradients

	for (unsigned int i = 0; i < lcvec.size(); i++) {
		lcvec[i].setElementPointer(&evec[i]);
		lcvec[i].setMaterialPointer(&mvec[i]);
		lcvec[i].setGlobalNodes(evec[i].getGlobalNodes());
	}

	std::vector<LocalConstructor> *lcvecptr;
	lcvecptr = &lcvec;

	Rod r1; // the rod object is constructed and fed all the vector pointers
	r1.setNvecP(nvecptr);
	r1.setFvecP(fvecptr);
	r1.setEvecP(evecptr);
	r1.setMvecP(mvecptr);
	r1.setLCvecP(lcvecptr);
	r1.setParallel(false); // not solving in parallel
	r1.Update(2);
	r1.setrho(0.0); //setting initial gravity to zero
	r1.calcGlobaltoFree();

	std::vector<Inextensibility> Inext(fvecptr->size()); // constructing the inextensibility constraint

	for (unsigned int i = 0; i < Inext.size(); i++) {
		Inext[i].setFramePointer(&(fvec[i])); //inextensibility is a frame constraint
		Inext[i].setdS(si[i+1] - si[i]);
		Inext[i].setGlobalNodes(fvec[i].getGlobalNodes());
	}
	std::vector<Constraint*> ConstP;
	for (unsigned int i = 1; i < (Inext.size()); i++) {
		ConstP.push_back(&Inext[i]); // adding inextensibility to the constraint pointer vector

	}

	std::vector<Constraint*> *PConstP;
	PConstP = &ConstP; // pointer to the vector of pointer constraints


	r1.setPConstP(PConstP);
	r1.Update(2); // Hessian level update to initialize everything.

	Solver s1; // initializing the solver object
	s1.setRodP(&r1); //feeding in a pointer to the rod object
	s1.setSolveDetailLevel(SolverRead);

	std::vector<double> mult1 = Linspace(0.0, 1.0, Nstep); // setting up solve iterations
	cout << "Bending Modulus: " << BMod[0] << ", " << BMod[1] << ", " << BMod[2] << endl;


	// Solve path. I.e. solving for progressive solutions
	for (unsigned int i = 0; i < Nstep; i++) {
		cout << "Step: " << i << endl;

		Eigen::Vector3d k0temp;
		k0temp.setZero();
		k0temp[0] = mult1[i]*1.0*pi; // setting curvature vector
		for (unsigned int j = 0; j < mvec.size(); j++){
			mvec[j].setkappa0(k0temp); // feeding in the curvature vector to the material class
		}

		double rhotemp = mult1[i] * rho; // setting global gravity parameter
		Eigen::Vector3d fendtemp;
		fendtemp = mult1[i]*fend;
		nvec.back().setFnode(fendtemp);


		r1.setrho((si[1] - si[0])*rhotemp); // Updating the global gravity parameter
		r1.Update(2); // Hessian level update
		s1.Solve(); // Newton-Raphson iteration to solve for equilibrium solution
		r1.setdtoD(); // Update d from D


		std::string directoryname;
		std::ostringstream ss;
		ss << foldername << "/" << filename << std::setw(5) << std::setfill('0')
						<< i << ".vtk";

		writevtk(ss.str(), nvec, fvec); // Saving configuration as a .vtk file.
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
		myfile << nvec[i].getX()[0]  << " "
				<< nvec[i].getX()[1] << " "
				<< nvec[i].getX()[2]  << "\n";
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

