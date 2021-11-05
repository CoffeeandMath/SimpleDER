/*
 * X3Const.h
 *
 *  Created on: Jul 9, 2020
 *      Author: KevinKorner
 */

#ifndef X3CONST_H_
#define X3CONST_H_
#include<Constraint.h>
#include<Node.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

class X3Const: public Constraint {
public:
	X3Const();
	~X3Const() {
	}
	;
	void setNodePointer(Node*);
	void Update(int);
	void setX3min(double);

private:
	void calcC();
	void calcdC();
	void calcddC();
	Node *NodeP;
	double X3min;

};

#endif /* SRC_X3CONST_H_ */
