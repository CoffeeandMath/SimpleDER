/*
 * Inextensibility.h
 *
 *  Created on: Jun 19, 2020
 *      Author: KevinKorner
 */

#ifndef INEXTENSIBILITY_H_
#define INEXTENSIBILITY_H_
#include<Constraint.h>
#include<Frame.h>


class Inextensibility: public Constraint {
public:
	Inextensibility();
	~Inextensibility(){};
	void Update(int);
	void setFramePointer(Frame*);
	void setdS(double);
protected:

private:
	void calcC();
	void calcdC();
	void calcddC();

	Frame* FrameP;

	double dS;
	double dSinv;
};

#endif /* INEXTENSIBILITY_H_ */
