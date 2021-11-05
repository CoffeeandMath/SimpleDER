/*
 * X3Const.cpp
 *
 *  Created on: Jul 9, 2020
 *      Author: KevinKorner
 */

#include "X3Const.h"

X3Const::X3Const() {
	NodeP = NULL;
	dC.resize(4);
	dC.setZero();
	ddC.resize(4, 4);
	ddC.setZero();
	X3min = -10000.0;
}

void X3Const::Update(int Order) {
	calcC();
	if (Order > 0) {
		calcdC();
	}
	if (Order > 1) {
		calcddC();
	}
}

void X3Const::calcC() {
	C = NodeP->getX()[2] - X3min;
}
void X3Const::calcdC() {

	dC[2] = 1.0;

}
void X3Const::calcddC() {

}
void X3Const::setNodePointer(Node *np) {
	NodeP = np;
}
void X3Const::setX3min(double xmin) {
	X3min = xmin;
}
