/*
 * Rod.h
 *
 *  Created on: Jun 16, 2020
 *      Author: kevin
 */

#ifndef ROD_H_
#define ROD_H_

#include<Node.h>
#include<Frame.h>
#include<Element.h>
#include<Material.h>
#include<LocalConstructor.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Constraint.h>
#include "omp.h"
#include <algorithm>
#include <iterator>
#include <iostream>
using namespace std;
typedef Eigen::Triplet<double> T;
class Rod {
public:
	Rod();

	double getE();
	double getEbend();
	double getEext();
	double getEgrav();
	Eigen::VectorXd getdE();
	Eigen::VectorXd getdEbend();
	Eigen::VectorXd getdEext();
	Eigen::VectorXd getdEgrav();

	std::vector<T> getddE();
	Eigen::SparseMatrix<double> getddEbend();
	void setrho(double);
	double getrho();
	void setNn(int);
	int getNn();

	void Update(int);

	void setNvecP(std::vector<Node>*);
	void setFvecP(std::vector<Frame>*);
	void setEvecP(std::vector<Element>*);
	void setMvecP(std::vector<Material>*);
	void setLCvecP(std::vector<LocalConstructor>*);
	void setParallel(bool);

	void setGravDir(Eigen::Vector3d&);

	void calcGlobaltoFree();

	std::vector<int> getGlobaltoFree();
	std::vector<int> getFreetoGlobal();

	Eigen::VectorXd ExtractGenCoord();
	Eigen::VectorXd ExtractGenCoordAll();

	void ApplyGenCoordAll(Eigen::VectorXd&);
	void ApplyGenCoord(Eigen::VectorXd&);
	void setdtoD();
	void setPConstP(std::vector<Constraint*>*);
	std::vector<Constraint*>* getPConstP();

	int getNc();

protected:

private:
	void calcE();
	void calcdE();
	void calcddE();
	void calcEbend();
	void calcEext();
	void calcEgrav();

	void calcdEbend();
	void calcdEext();
	void calcdEgrav();

	void calcddEbend();

	void calcNdof();


	int Nn;
	int Nc;
	int Ndof;

	std::vector<Node>*NvecP;
	std::vector<Frame>*FvecP;
	std::vector<Element>*EvecP;
	std::vector<Material>*MvecP;
	std::vector<LocalConstructor>*LCvecP;

	std::vector<Constraint*>* PConstP;

	std::vector<int> GlobaltoFree;
	std::vector<int> FreetoGlobal;

	double E;
	double Ebend;
	double Eext;
	double Egrav;

	Eigen::VectorXd dE;
	Eigen::VectorXd dEbend;
	Eigen::VectorXd dEext;
	Eigen::VectorXd dEgrav;

	Eigen::SparseMatrix<double> ddE;
	Eigen::SparseMatrix<double> ddEbend;

	double rho;
	Eigen::Vector3d gravdir;

	bool Parallelized;

	Eigen::VectorXd GenCoord;
	Eigen::VectorXd GenCoordAll;

	void calcGenCoord();

	std::vector<T> ddEbendtrips;
	std::vector<T> ddEtrips;
	template < typename T>
	std::pair<bool, int > findInVector(const std::vector<T>  & , const T  & );

};

#endif /* ROD_H_ */
