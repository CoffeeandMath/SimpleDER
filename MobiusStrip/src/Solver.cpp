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
	mu = 1000.0;
	Nactive = 0;
	firstsolve = true;
	fk = 0.0;
	Nfree = 0;
	N = 0;
	SolveDetail = 0;
	kiters = 0;
	tol = 10;
	tsolve = 0.0;
	laststep = 0.0;
	LastInRemoved = -1;
}

void Solver::Solve() {
	//ALSolve();

	if (RodP->getNcIn() == 0) {
		SQPSolve();
	} else {
		SQPASSolve();
	}

}

void Solver::SQPSolve() {
	double eta = 0.0;
	double tau = 0.5;
	double sigma = 1.0;
	double rho = 0.5;
	RodP->Update(2);
	Eigen::VectorXd xk = RodP->ExtractGenCoord();
	Eigen::VectorXd xkall = RodP->ExtractGenCoordAll();
	Nfree = xk.size();

	N = xkall.size();
	mu = 1.0 * N;
	unsigned int k = 0;
	tol = 10;
	double tol0 = pow(10.0, -10.0) * N;
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

		solver.compute(M);
		if (solver.info() != Eigen::Success) {
			cout << "decomposition failed" << endl;
			return;
		}
		Eigen::VectorXd p = solver.solve(v);
		if (solver.info() != Eigen::Success) {
			cout << "solving failed" << endl;
			return;
		}
		double alpha = 1.0;
		for (unsigned int i = 0; i < xk.size(); i++) {
			xk[i] += alpha * p[i];
		}
		for (unsigned int i = 0; i < RodP->getNc(); i++) {
			lambda[i] = alpha * p[Nfree + i];
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
	auto duration = duration_cast < microseconds > (stop - start);
	tsolve = duration.count();
	SolveSQPReadout();

}

void Solver::SQPASSolve() {
	double tau = 0.5;
	RodP->Update(2);
	Eigen::VectorXd xk = RodP->ExtractGenCoord();
	Eigen::VectorXd xkall = RodP->ExtractGenCoordAll();
	Nfree = xk.size();

	N = xkall.size();
	mu = 1.0 * N;
	unsigned int k = 0;
	tol = 10;
	double tol0 = pow(10.0, -8.0) * N;
	unsigned int kmax = 30;

	if (lambda.size() != RodP->getNc()) {
		lambda.resize(RodP->getNc());
		lambda.setZero();

	}

	if (isActive.size() != RodP->getNcIn()) {
		cout << "resizing LambdaIn" << endl;
		isActive.resize(RodP->getNcIn());
		lambdaIn.resize(RodP->getNcIn());
		for (unsigned int i = 0; i < RodP->getNcIn(); i++) {
			isActive[i] = false;
		}
		lambdaIn.setZero();
	}
	Nactive = 0;
	for (unsigned int i = 0; i < RodP->getNcIn(); i++) {
		if (isActive[i]) {
			Nactive++;
		}
	}

	auto start = high_resolution_clock::now();
	if (Nactive > 0) {
		UpdateSQPAStrips();
	} else {
		UpdateSQPtrips();
	}

	while ((tol > tol0) && (k < kmax)) {
		k++;

		Nactive = 0;

		for (unsigned int i = 0; i < RodP->getNcIn(); i++) {
			if (isActive[i]) {
				Nactive++;
			}
		}
		SpMat M(Nfree + RodP->getNc() + Nactive,
				Nfree + RodP->getNc() + Nactive);
		Eigen::VectorXd v(Nfree + RodP->getNc() + Nactive);

		v << (-1.0) * Dfk, ck;
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
		/*
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
		 */
		//solver.analyzePattern(M);
		solver.compute(M);
		if (solver.info() != Eigen::Success) {
			cout << "decomposition failed" << endl;
			return;
		}
		Eigen::VectorXd p(v.size());
		p = solver.solve(v);

		if (solver.info() != Eigen::Success) {
			cout << "solving failed" << endl;
			return;
		}

		std::vector<int> activeindices;
		for (unsigned int i = 0; i < RodP->getNcIn(); i++) {
			if (isActive[i]) {
				activeindices.push_back(i);
			}
		}
		RodP->SetActiveIndices(activeindices);

		Eigen::VectorXd px(Nfree);
		for (unsigned int i = 0; i < xk.size(); i++) {
			px[i] = p[i];
		}
		for (unsigned int i = 0; i < RodP->getNc(); i++) {
			lambda[i] = p[Nfree + i];
		}
		for (unsigned int i = 0; i < activeindices.size(); i++) {
			lambdaIn[activeindices[i]] = p[Nfree + RodP->getNc() + i];
		}

		Eigen::VectorXd xktemp1(p.size());
		double alphakls = 1.0;
		double tauls = 0.5;

		xktemp1 = xk + px;
		laststep = px.norm();
		RodP->ApplyGenCoord(xktemp1);
		RodP->Update(0);
		bool NanDetected = false;
		while (isnan(RodP->getE())) {
			alphakls *= tauls;
			xktemp1 = xk + alphakls * px;
			RodP->ApplyGenCoord(xktemp1);
			RodP->Update(0);
			NanDetected = true;
		}
		if (NanDetected) {
			p *= alphakls;
		}

		double alphak = 1.0;
		tol = px.norm();
		if (tol < tol0) {
			int cnttemp = 0;

			/*for (unsigned int i = (Nfree + RodP->getNc()); i < p.size(); i++) {
			 while (~isActive[cnttemp]) {
			 cnttemp++;
			 }
			 lambdaIn[cnttemp] = lambda[i];
			 cnttemp++;
			 }
			 */

			int smallestNegLambdaIndex = -1;
			double lambdaNegSmall = 0.0;
			for (unsigned int i = 0; i < lambdaIn.size(); i++) {
				if (isActive[i] && (lambdaIn[i] < 0.0)) {
					if (lambdaIn[i] < lambdaNegSmall) {
						smallestNegLambdaIndex = i;
						lambdaNegSmall = lambdaIn[i];
					}
				}
			}
			for (unsigned int i = 0; i < xk.size(); i++) {
				xk[i] += px[i];
			}
			tol = px.norm();
			laststep = tol;
			if (smallestNegLambdaIndex > -1) {
				isActive[smallestNegLambdaIndex] = false;
				LastInRemoved = smallestNegLambdaIndex;
				lambdaIn[smallestNegLambdaIndex] = 0.0;
				tol = 10.0;
				//kmax += 10;
				cout << "Removed constraint: " << smallestNegLambdaIndex
						<< endl;
				//k = 0;
			}
		} else {

			int blockingindex = -1;
			//Compute alphak
			Eigen::VectorXd xktemp;
			xktemp = xk + alphak * px;
			RodP->ApplyGenCoord(xktemp);
			RodP->Update(0);
			for (unsigned int i = 0; i < isActive.size(); i++) {
				if (!isActive[i] && i != LastInRemoved) {
					int itermax = 0;
					if ((((*(RodP->getPConstInP()))[i])->getC()) < 0.0) {
						blockingindex = i;
						while (((((*(RodP->getPConstInP()))[i])->getC())
								< -0.0001) && itermax < 15) {

							itermax++;
							alphak *= tau;
							xktemp = xk + alphak * px;
							RodP->ApplyGenCoord(xktemp);
							RodP->Update(0);

						}
					}
				}
			}

			if (blockingindex > -1) {
				//k = 0;
				isActive[blockingindex] = true;
				//kmax += 10;
				cout << "Added constraint: " << blockingindex << endl;
				tol = 10.0;
			}

		}
		Eigen::VectorXd xktemp;
		xk += alphak * px;
		RodP->ApplyGenCoord(xk);
		RodP->Update(2);
		laststep = alphak * px.norm();
		//cout << tol / tol0 << endl;
		Nactive = 0;
		for (unsigned int i = 0; i < RodP->getNcIn(); i++) {
			if (isActive[i]) {
				Nactive++;
			}
		}
		//Calculate update
		if (Nactive > 0) {
			UpdateSQPAStrips();
		} else {
			UpdateSQPtrips();
		}
	}

	kiters = k;
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast < microseconds > (stop - start);
	tsolve = duration.count();
	SolveSQPReadout();
}

void Solver::ALSolve() {
	RodP->Update(2);
	Eigen::VectorXd xk = RodP->ExtractGenCoord();
	Eigen::VectorXd xkall = RodP->ExtractGenCoordAll();
	Nfree = xk.size();

	N = xkall.size();
	if (mu < 10.0) {
		mu = 10.0 * N;
	}
	unsigned int k = 0;
	tol = 10;
	double tol0 = pow(10.0, -12.0) * N;
	unsigned int kmax = 30;

	if (lambda.size() != RodP->getNc()) {
		lambda.resize(RodP->getNc());
		lambda.setZero();

	}
	auto start = high_resolution_clock::now();

	while ((tol > tol0) && (k < kmax)) {
		k++;
		//cout << "AL Iterate" << endl;

		NewtonRaphsonSolve();

		tol = c1norm();
		cout << "tol value:" << tol << endl;
		for (unsigned int i = 0; i < RodP->getNc(); i++) {
			lambda[i] -= mu * (((*(RodP->getPConstP()))[i])->getC());
		}
		//mu = 10.0 * mu;

	}

	kiters = k;
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast < microseconds > (stop - start);
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

	for (unsigned int i = 0; i < RodP->getNc(); i++) {
		//Dfk[(RodP->getGlobaltoFree())[]]
		for (unsigned int j = 0;
				j < ((*(RodP->getPConstP()))[i]->getGlobalNodes()).size();
				j++) {
			if (RodP->getGlobaltoFree()[(*(RodP->getPConstP()))[i]->getGlobalNodes()[j]]
					>= 0) {
				Dfk[RodP->getGlobaltoFree()[(*(RodP->getPConstP()))[i]->getGlobalNodes()[j]]] +=
						mu * (*(RodP->getPConstP()))[i]->getC()
								* ((*(RodP->getPConstP()))[i]->getdC()[j]);
			}

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
			if (RodP->getGlobaltoFree()[(*(RodP->getPConstP()))[i]->getGlobalNodes()[j]]
					>= 0) {
				Aktrips.push_back(
						T(i,
								RodP->getGlobaltoFree()[(*(RodP->getPConstP()))[i]->getGlobalNodes()[j]],
								(*(RodP->getPConstP()))[i]->getdC()[j]));
			}

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

		Mtrips.push_back(
				T(Aktrips[i].row() + Nfree, Aktrips[i].col(),
						-1.0 * Aktrips[i].value()));
		Mtrips.push_back(
				T(Aktrips[i].col(), Aktrips[i].row() + Nfree,
						-1.0 * Aktrips[i].value()));

	}
}

void Solver::UpdateSQPAStrips() {

	fk = RodP->getE();

	std::vector<T> DDLtemp = RodP->getddE();

	std::vector<int> activeindices;
	for (unsigned int i = 0; i < RodP->getNcIn(); i++) {
		if (isActive[i]) {
			activeindices.push_back(i);
		}
	}
	RodP->SetActiveIndices(activeindices);
	ck.resize(RodP->getNc() + activeindices.size());
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

	for (unsigned int i = 0; i < activeindices.size(); i++) {

		int mnsize =
				((*(RodP->getPConstInP()))[activeindices[i]])->getGlobalNodes().size();
		for (unsigned int m = 0; m < mnsize; m++) {
			for (unsigned int n = 0; n < mnsize; n++) {

				DDLtemp.push_back(
						T(
								(((*(RodP->getPConstInP()))[activeindices[i]])->getGlobalNodes())[m],
								(((*(RodP->getPConstInP()))[activeindices[i]])->getGlobalNodes())[n],
								-1.0 * lambdaIn(activeindices[i])
										* (((*(RodP->getPConstInP()))[activeindices[i]])->getddC())(
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

	for (unsigned int i = 0; i < activeindices.size(); i++) {

		int mnsize =
				((*(RodP->getPConstInP()))[activeindices[i]])->getGlobalNodes().size();
		for (unsigned int m = 0; m < mnsize; m++) {
			for (unsigned int n = 0; n < mnsize; n++) {

				DDLtemp.push_back(
						T(
								(((*(RodP->getPConstInP()))[activeindices[i]])->getGlobalNodes())[m],
								(((*(RodP->getPConstInP()))[activeindices[i]])->getGlobalNodes())[n],
								mu
										* ((((*(RodP->getPConstInP()))[activeindices[i]])->getC())
												* ((((*(RodP->getPConstInP()))[activeindices[i]])->getddC())(
														m, n))
												+ ((((*(RodP->getPConstInP()))[activeindices[i]])->getdC())[m])
														* ((((*(RodP->getPConstInP()))[activeindices[i]])->getdC())[n]))));
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

	for (unsigned int i = 0; i < RodP->getNc(); i++) {
		//Dfk[(RodP->getGlobaltoFree())[]]
		for (unsigned int j = 0;
				j < ((*(RodP->getPConstP()))[i]->getGlobalNodes()).size();
				j++) {
			if (RodP->getGlobaltoFree()[(*(RodP->getPConstP()))[i]->getGlobalNodes()[j]]
					>= 0) {
				Dfk[RodP->getGlobaltoFree()[(*(RodP->getPConstP()))[i]->getGlobalNodes()[j]]] +=
						mu * (*(RodP->getPConstP()))[i]->getC()
								* ((*(RodP->getPConstP()))[i]->getdC()[j]);
			}

		}
	}

	for (unsigned int i = 0; i < activeindices.size(); i++) {
		//Dfk[(RodP->getGlobaltoFree())[]]

		for (unsigned int j = 0;
				j
						< ((*(RodP->getPConstInP()))[activeindices[i]]->getGlobalNodes()).size();
				j++) {
			if (RodP->getGlobaltoFree()[(*(RodP->getPConstInP()))[activeindices[i]]->getGlobalNodes()[j]]
					>= 0) {
				Dfk[RodP->getGlobaltoFree()[(*(RodP->getPConstInP()))[activeindices[i]]->getGlobalNodes()[j]]] +=
						mu * (*(RodP->getPConstInP()))[activeindices[i]]->getC()
								* ((*(RodP->getPConstInP()))[activeindices[i]]->getdC()[j]);
			}

		}
	}

	Aktrips.clear();
	if (Aktrips.capacity() != 20 * (RodP->getNc()) + Nactive) {
		Aktrips.reserve(20 * (RodP->getNc() + Nactive));
	}

	if (ck.size() != (RodP->getNc() + Nactive)) {
		ck.resize(RodP->getNc() + Nactive);
	}

	for (unsigned int i = 0; i < (RodP->getNc()); i++) {
		ck[i] = (*(RodP->getPConstP()))[i]->getC();
		for (unsigned int j = 0;
				j < ((*(RodP->getPConstP()))[i]->getGlobalNodes()).size();
				j++) {
			if (RodP->getGlobaltoFree()[(*(RodP->getPConstP()))[i]->getGlobalNodes()[j]]
					>= 0) {
				Aktrips.push_back(
						T(i,
								RodP->getGlobaltoFree()[(*(RodP->getPConstP()))[i]->getGlobalNodes()[j]],
								(*(RodP->getPConstP()))[i]->getdC()[j]));
			}

		}

	}

	for (unsigned int i = 0; i < activeindices.size(); i++) {
		ck[RodP->getNc() + i] =
				(*(RodP->getPConstInP()))[activeindices[i]]->getC();
		for (unsigned int j = 0;
				j
						< ((*(RodP->getPConstInP()))[activeindices[i]]->getGlobalNodes()).size();
				j++) {
			if (RodP->getGlobaltoFree()[(*(RodP->getPConstInP()))[activeindices[i]]->getGlobalNodes()[j]]
					>= 0) {
				Aktrips.push_back(
						T(RodP->getNc() + i,
								RodP->getGlobaltoFree()[(*(RodP->getPConstInP()))[activeindices[i]]->getGlobalNodes()[j]],
								(*(RodP->getPConstInP()))[activeindices[i]]->getdC()[j]));
			}

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

		Mtrips.push_back(
				T(Aktrips[i].row() + Nfree, Aktrips[i].col(),
						-1.0 * Aktrips[i].value()));
		Mtrips.push_back(
				T(Aktrips[i].col(), Aktrips[i].row() + Nfree,
						-1.0 * Aktrips[i].value()));

	}
}

double Solver::GetALval() {
	fk = RodP->getE();
	for (unsigned int i = 0; i < RodP->getNc(); i++) {
		fk += -lambda(i) * ((*(RodP->getPConstP()))[i]->getC())
				+ 0.5 * mu * pow((*(RodP->getPConstP()))[i]->getC(), 2.0);
	}
	return fk;
}

void Solver::UpdateALtrips() {
	fk = RodP->getE();

	Eigen::VectorXd Dfktemp = RodP->getdE();
	if (Dfk.size() != Nfree) {
		Dfk.resize(Nfree);
	}

	for (unsigned int i = 0; i < Dfktemp.size(); i++) {
		if ((RodP->getGlobaltoFree())[i] >= 0) {
			Dfk[(RodP->getGlobaltoFree())[i]] = Dfktemp[i];
		}
	}

	for (unsigned int i = 0; i < RodP->getNc(); i++) {
		//Dfk[(RodP->getGlobaltoFree())[]]
		for (unsigned int j = 0;
				j < ((*(RodP->getPConstP()))[i]->getGlobalNodes()).size();
				j++) {
			if (RodP->getGlobaltoFree()[(*(RodP->getPConstP()))[i]->getGlobalNodes()[j]]
					>= 0) {
				Dfk[RodP->getGlobaltoFree()[(*(RodP->getPConstP()))[i]->getGlobalNodes()[j]]] +=
						mu * (*(RodP->getPConstP()))[i]->getC()
								* ((*(RodP->getPConstP()))[i]->getdC()[j]);
				Dfk[RodP->getGlobaltoFree()[(*(RodP->getPConstP()))[i]->getGlobalNodes()[j]]] -=
						lambda(i) * ((*(RodP->getPConstP()))[i]->getdC()[j]);
			}

		}
	}

}

void Solver::setSolveDetailLevel(int in) {
	SolveDetail = in;
}

void Solver::SolveSQPReadout() {
	if (SolveDetail == 0) {

	} else if (SolveDetail == 1) {
		cout << "Last Step Size: " << laststep << endl;
		cout << "Number of Iterations: " << kiters << endl;
		cout << "Time to solve (seconds): "
				<< ((double) tsolve) / pow(10.0, 6.0) << endl;
		cout << "Number of Active Constraints: " << Nactive << endl;

	}

}

double Solver::c1norm() {
	double ctemp = 0.0;
	for (unsigned int i = 0; i < RodP->getNc(); i++) {
		ctemp += abs((*(RodP->getPConstP()))[i]->getC());
	}
	return ctemp;
}

void Solver::NewtonRaphsonSolve() {
	Eigen::VectorXd xk = RodP->ExtractGenCoord();
	Eigen::VectorXd p;

	UpdateALtrips();
	double dfnorm = Dfk.norm();
	double cnorm = c1norm();
	int k = 0;
	int kmax = 150;
	double tolNR = pow(10.0, -8.0) * Dfk.size();
	while (k < kmax && dfnorm > tolNR) {
		//cout << "solving" << endl;
		double alpha = 1.0;
		k++;

		Eigen::VectorXd v(Dfk.size());
		v << (-1.0) * Dfk;
		/*
		 SpMat Bk(Nfree, Nfree);

		 Bk.setFromTriplets(Bktrips.begin(), Bktrips.end());
		 Bk.makeCompressed();

		 solver.compute(Bk);
		 if (solver.info() != Eigen::Success) {
		 cout << "decomposition failed" << endl;
		 return;
		 }
		 p = solver.solve(v);
		 if (solver.info() != Eigen::Success) {
		 cout << "solving failed" << endl;
		 return;
		 }
		 */
		double Eold = GetALval();
		p = alpha * v;
		Eigen::VectorXd xktemp;
		xktemp = xk + p;
		RodP->ApplyGenCoord(xktemp);
		RodP->Update(0);
		double Enew = GetALval();
		//p = 0.000001 * v;

		while (Enew > (Eold) && v.norm() > pow(10.0, -16.0)) {
			alpha *= 0.5;
			p = alpha * v;
			//cout << p.norm() << endl;
			xktemp = xk + p;
			RodP->ApplyGenCoord(xktemp);
			RodP->Update(0);
			Enew = GetALval();
			//cout << "Enew: " << Enew << ", Eold: " << Eold << endl;
		}

		xk += 0.5 * p;
		RodP->ApplyGenCoord(xk);
		RodP->Update(2);
		UpdateALtrips();
		dfnorm = p.norm();
		//cout << dfnorm << endl;
	}
	cout << k << endl;

}
