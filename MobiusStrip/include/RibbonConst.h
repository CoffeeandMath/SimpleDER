/*
 * RibbonConst.h
 *
 *  Created on: Jul 1, 2020
 *      Author: KevinKorner
 */

#ifndef RIBBONCONST_H_
#define RIBBONCONST_H_
#include<Constraint.h>
#include<Element.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
class RibbonConst: public Constraint {
public:
	RibbonConst();
	~RibbonConst() {
	}
	;
	void Update(int);
	void setElementPointer(Element*);
	void setdS(double);

private:
	void calcC();
	void calcdC();
	void calcddC();

	Element *ElemP = NULL;

	double dS;
	double dSinv;
};

#endif /* RIBBONCONST_H_ */
