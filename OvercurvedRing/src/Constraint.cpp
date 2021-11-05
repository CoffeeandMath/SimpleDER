/*
 * Constraint.cpp
 *
 *  Created on: Jun 19, 2020
 *      Author: KevinKorner
 */

#include "Constraint.h"

Constraint::Constraint() {
	// TODO Auto-generated constructor stub
	C = 0;

}
double Constraint::getC()
{
	return C;
}
Eigen::VectorXd Constraint::getdC()
{
	return dC;
}
Eigen::MatrixXd Constraint::getddC()
{
	return ddC;
}
void Constraint::setGlobalNodes(Eigen::VectorXi GN)
{
    globalNodes = GN;
}
Eigen::VectorXi Constraint::getGlobalNodes()
{
    return globalNodes;
}
