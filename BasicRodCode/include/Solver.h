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
	Rod *RodP;
	Eigen::VectorXd lambda;
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

	void UpdateSQPtrips();
	int Nfree;
	int N;

	int SolveDetail;
	void SolveSQPReadout();

	int kiters;

	double tol;
	double tsolve;

};

#endif /* SOLVER_H_ */
