/*
 * EtapLimP.h
 *
 *  Created on: Jul 9, 2020
 *      Author: KevinKorner
 */

#ifndef ETAPLIMP_H_
#define ETAPLIMP_H_
#include<Constraint.h>
#include<Material.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

class EtapLimP: public Constraint {
public:
	EtapLimP();
	~EtapLimP() {
	}
	;
	void setMaterialPointer(Material*);
	void Update(int);
	void setBarrierLoc(double);


private:
	void calcC();
	void calcdC();
	void calcddC();
	Material *MatP;
	double Bloc;


};

#endif /* SRC_X3CONST_H_ */
