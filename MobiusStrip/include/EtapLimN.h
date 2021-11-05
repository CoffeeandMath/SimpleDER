/*
 * EtapLimP.h
 *
 *  Created on: Jul 9, 2020
 *      Author: KevinKorner
 */

#ifndef ETAPLIMN_H_
#define ETAPLIMN_H_
#include<Constraint.h>
#include<Material.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

class EtapLimN: public Constraint {
public:
	EtapLimN();
	~EtapLimN() {
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
