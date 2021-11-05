/*
 * Solver.cpp
 *
 *  Created on: Jun 22, 2020
 *      Author: KevinKorner
 */

#include "Solver.h"

Solver::Solver() {
	// TODO Auto-generated constructor stub

	RodP = NULL;
	mu = 0;
	Eigen::setNbThreads(30);
	Eigen::initParallel();
	firstsolve = true;
	fk = 0.0;
	Nfree = 0;
	N = 0;
	SolveDetail = 0;
	kiters = 0;
	tol = 10;
	tsolve = 0.0;
}

void Solver::Solve() {

	SQPSolve();

}

void Solver::SQPSolve() {
	double eta = 0.0;
	double tau = 0.5;
	double sigma = 1.0;
	double rho = 0.5;

	Eigen::VectorXd xk = RodP->ExtractGenCoord();
	Eigen::VectorXd xkall = RodP->ExtractGenCoordAll();
	Nfree = xk.size();

	N = xkall.size();
	mu = 1.0 * N;
	unsigned int k = 0;
	tol = 10;
	double tol0 = pow(10.0, -12.0) * N;
	unsigned int kmax = 15;

	if (lambda.size() != RodP->getNc()) {
		lambda.resize(RodP->getNc());
		lambda.setZero();

	}
	auto start = high_resolution_clock::now();

	UpdateSQPtrips();

	while ((tol > tol0) && (k < kmax)) {
		k++;
		Eigen::VectorXd v(Dfk.size() + ck.size());
		v << (-1.0) * Dfk, ck;

		SpMat M(Nfree + RodP->getNc(), Nfree + RodP->getNc());
		/*
		 cout << Mtrips.size() << endl;
		 for (unsigned int i = 0; i < Mtrips.size(); i++)
		 {
		 cout << "i = " << i << ": " << Mtrips[i].row() << ",  " << Mtrips[i].col() << ", " << Mtrips[i].value() << endl;
		 }
		 */

		M.setFromTriplets(Mtrips.begin(), Mtrips.end());
		M.makeCompressed();

		//cout << "got here" << endl;

		// solve Ax = b
		if (firstsolve) {
			cout << "Analyzing Hessian Structure" << endl;
			solver.analyzePattern(M);
			firstsolve = false;
		}

		solver.factorize(M);
		if (solver.info() != Eigen::Success) {
			cout << "decomposition failed" << endl;
			return;
		}
		Eigen::VectorXd p = solver.solve(v);
		if (solver.info() != Eigen::Success) {
			cout << "solving failed" << endl;
			return;
		}

		for (unsigned int i = 0; i < xk.size(); i++) {
			xk[i] += p[i];
		}
		for (unsigned int i = 0; i < RodP->getNc(); i++) {
			lambda[i] = p[Nfree + i];
		}

		RodP->ApplyGenCoord(xk);
		RodP->Update(2);
		tol = 0;
		for (unsigned int i = 0; i < Nfree; i++) {
			tol += pow(p[i], 2.0);
		}

		//cout << tol / tol0 << endl;

		//Calculate update
		UpdateSQPtrips();

	}

	kiters = k;
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	tsolve = duration.count();
	SolveSQPReadout();

}

void Solver::setRodP(Rod *RP) {
	RodP = RP;
}

Eigen::MatrixXd Solver::TripList2EMat(std::vector<T> tlist) {
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

void Solver::UpdateSQPtrips() {
	fk = RodP->getE();
	std::vector<T> DDLtemp = RodP->getddE();

	for (unsigned int i = 0; i < RodP->getNc(); i++) {
		int mnsize = ((*(RodP->getPConstP()))[i])->getGlobalNodes().size();
		for (unsigned int m = 0; m < mnsize; m++) {
			for (unsigned int n = 0; n < mnsize; n++) {

				DDLtemp.push_back(
						T((((*(RodP->getPConstP()))[i])->getGlobalNodes())[m],
								(((*(RodP->getPConstP()))[i])->getGlobalNodes())[n],
								-1.0 * lambda(i)
										* (((*(RodP->getPConstP()))[i])->getddC())(
												m, n)));
			}

		}
	}

	for (unsigned int i = 0; i < RodP->getNc(); i++) {
		int mnsize = ((*(RodP->getPConstP()))[i])->getGlobalNodes().size();
		for (unsigned int m = 0; m < mnsize; m++) {
			for (unsigned int n = 0; n < mnsize; n++) {

				DDLtemp.push_back(
						T((((*(RodP->getPConstP()))[i])->getGlobalNodes())[m],
								(((*(RodP->getPConstP()))[i])->getGlobalNodes())[n],
								mu
										* ((((*(RodP->getPConstP()))[i])->getC())
												* ((((*(RodP->getPConstP()))[i])->getddC())(
														m, n))
												+ ((((*(RodP->getPConstP()))[i])->getdC())[m])
														* ((((*(RodP->getPConstP()))[i])->getdC())[n]))));
			}

		}
	}
	Bktrips.clear();
	if (Bktrips.capacity() != DDLtemp.size()) {
		Bktrips.reserve(DDLtemp.size());
	}
	for (unsigned int i = 0; i < DDLtemp.size(); i++) {
		if ((RodP->getGlobaltoFree())[DDLtemp[i].row()] >= 0
				&& (RodP->getGlobaltoFree())[DDLtemp[i].col()] >= 0) {
			Bktrips.push_back(
					T((RodP->getGlobaltoFree())[DDLtemp[i].row()],
							(RodP->getGlobaltoFree())[DDLtemp[i].col()],
							DDLtemp[i].value()));
		}
	}

	Eigen::VectorXd Dfktemp = RodP->getdE();
	if (Dfk.size() != Nfree) {
		Dfk.resize(Nfree);
	}

	for (unsigned int i = 0; i < Dfktemp.size(); i++) {
		if ((RodP->getGlobaltoFree())[i] >= 0) {
			Dfk[(RodP->getGlobaltoFree())[i]] = Dfktemp[i];
		}
	}

	Aktrips.clear();
	if (Aktrips.capacity() != 20 * (RodP->getNc())) {
		Aktrips.reserve(20 * (RodP->getNc()));
	}

	if (ck.size() != RodP->getNc()) {
		ck.resize(RodP->getNc());
	}

	for (unsigned int i = 0; i < (RodP->getNc()); i++) {
		ck[i] = (*(RodP->getPConstP()))[i]->getC();
		for (unsigned int j = 0;
				j < ((*(RodP->getPConstP()))[i]->getGlobalNodes()).size();
				j++) {
			Aktrips.push_back(
					T(i, (*(RodP->getPConstP()))[i]->getGlobalNodes()[j],
							(*(RodP->getPConstP()))[i]->getdC()[j]));
		}

	}

	Mtrips.clear();
	if (Mtrips.capacity() != (Bktrips.size() + 2 * Aktrips.size())) {
		Mtrips.reserve(Bktrips.size() + 2 * Aktrips.size());
	}

	for (unsigned int i = 0; i < Bktrips.size(); i++) {
		Mtrips.push_back(Bktrips[i]);
	}

	for (unsigned int i = 0; i < Aktrips.size(); i++) {

		if ((RodP->getGlobaltoFree())[Aktrips[i].col()] >= 0) {

			Mtrips.push_back(
					T(Aktrips[i].row() + Nfree,
							RodP->getGlobaltoFree()[Aktrips[i].col()],
							-1.0 * Aktrips[i].value()));
			Mtrips.push_back(
					T(RodP->getGlobaltoFree()[Aktrips[i].col()],
							Aktrips[i].row() + Nfree,
							-1.0 * Aktrips[i].value()));
		}
	}
}

void Solver::setSolveDetailLevel(int in) {
	SolveDetail = in;
}

void Solver::SolveSQPReadout() {
	if (SolveDetail == 0) {

	} else if (SolveDetail == 1) {
		cout << "Last Step Size: " << tol << endl;
		cout << "Number of Iterations: " << kiters << endl;
		cout << "Time to solve (seconds): "
				<< ((double) tsolve) / pow(10.0, 6.0) << endl;

	}

}
