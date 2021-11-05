/*
 * k20.h
 *
 *  Created on: Jun 30, 2020
 *      Author: KevinKorner
 */

#ifndef K20_H_
#define K20_H_
#include<Constraint.h>
#include<Element.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
class k20: public Constraint {
public:
	k20();
	~k20() {
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

#endif /* K20_H_ */
