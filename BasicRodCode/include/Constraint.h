/*
 * Constraint.h
 *
 *  Created on: Jun 19, 2020
 *      Author: KevinKorner
 */

#ifndef CONSTRAINT_H_
#define CONSTRAINT_H_

#include <Eigen/Dense>

class Constraint {
public:
  Constraint();
	virtual ~Constraint(){};
	virtual void Update(int){};

	double getC();
	Eigen::VectorXd getdC();
	Eigen::MatrixXd getddC();
	void setGlobalNodes(Eigen::VectorXi GN);
	Eigen::VectorXi getGlobalNodes();

protected:
	double C;
	Eigen::VectorXd dC;
	Eigen::MatrixXd ddC;

	Eigen::VectorXi globalNodes;

private:


};

#endif /* CONSTRAINT_H_ */
