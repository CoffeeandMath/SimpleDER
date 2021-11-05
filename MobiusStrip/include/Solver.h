/*
 * Solver.h
 *
 *  Created on: Jun 22, 2020
 *      Author: KevinKorner
 */

#ifndef SOLVER_H_
#define SOLVER_H_
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>
//#include <Eigen/PardisoSupport>
#include <Eigen/Core>
#include <Rod.h>
#include <Constraint.h>
#include <math.h>
#include <chrono>
typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> SpMat;
using namespace std::chrono;

class Solver {
public:
	Solver();
	void Solve();
	void setRodP(Rod*);
	void setSolveDetailLevel(int);

private:
	//std::string method;
	void SQPSolve();
	void SQPASSolve();

	void ALSolve();
	void NewtonRaphsonSolve();
	Rod *RodP;
	Eigen::VectorXd lambda;
	Eigen::VectorXd lambdaIn;
	double mu;
	Eigen::MatrixXd TripList2EMat(std::vector<T>);
	Eigen::SimplicialLDLT<SpMat> solver;
	bool firstsolve;

	std::vector<T> Bktrips;
	std::vector<T> Mtrips;
	double fk;
	Eigen::VectorXd Dfk;
	Eigen::VectorXd ck;
	std::vector<T> Aktrips;

	std::vector<bool> isActive;

	void UpdateSQPtrips();
	void UpdateSQPAStrips();
	void UpdateALtrips();
	int Nfree;
	int N;

	int SolveDetail;
	void SolveSQPReadout();
	double GetALval();
	int LastInRemoved;

	int kiters;
	int Nactive;

	double tol;
	double laststep;
	double tsolve;
	double c1norm();

};

#endif /* SOLVER_H_ */
