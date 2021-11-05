#include "Frame.h"

Frame::Frame()
{
    phi = 0;
    eps = 0;
    t = {0,0,0};
    T = {0,0,0};
    d.Set(1,0,0,0);
    D.Set(1,0,0,0);
    ptw.Set(1,0,0,0);
    ppar.Set(1,0,0,0);
    B.setZero();
    depsdx.setZero();
    DB = Eigen::Tensor<double, 3> (3,3,3);
    DB.setZero();
    ddepsdxdx = Eigen::Tensor<double,4> (3,3,2,2);
    ddepsdxdx.setZero();
    FreePhibool = true;


}

void Frame::Update(int Order)
{

	calct();
	calcT();
	calcptw();
	calcppar();
	d = ppar*ptw*D;
	d.Normalize();

	if (Order > 0){
		B = (0.5/eps)*((1.0/(1.0+T.dot(t)))*(T*(T.cross(t).transpose()) - T.cross(t)*(t.transpose() + T.transpose())) + skw(T));
		for (unsigned int i = 0; i < 3; i++)
		{
			depsdx(i,0) = -t[i];
			depsdx(i,1) = t[i];
		}

	}

	if (Order > 1)
	{
		double f2i = 1.0 + T.dot(t);
		double epsinv = 1.0/eps;
		Eigen::Tensor<double, 3> dBideps(3,3,3);
		for (unsigned int i = 0; i < 3; i++)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				for (unsigned int k = 0; k < 3; k++)
				{
					dBideps(i,j,k) = -epsinv*B(i,j)*t[k];
				}
			}
		}

		Eigen::Tensor<double, 3> dBidf2(3,3,3);
		Eigen::Matrix3d M;
		Eigen::Vector3d Ticti = T.cross(t);
		M = T*Ticti.transpose() - Ticti*(t.transpose() + T.transpose());
		Eigen::Vector3d v;
		v = -1.0*(0.5*pow(epsinv,2.0))*(1.0/pow(f2i,2.0))*(T - t*(f2i-1));

		for (unsigned int i = 0; i < 3; i++)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				for (unsigned int k = 0; k < 3; k++)
				{
					dBidf2(i,j,k) = M(i,j)*v[k];
				}
			}
		}

		Eigen::Tensor<double,3> dBidti(3,3,3);
		Eigen::Matrix3d M1;
		M1 = skw(T) - Ticti*t.transpose();
		Eigen::Matrix3d tiprojeps;
		tiprojeps = Eigen::Matrix3d::Identity() - t*t.transpose();

		double sc = (0.5*pow(epsinv,2.0)/f2i);
		for (unsigned int i = 0; i < 3; i++)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				for (unsigned int k = 0; k < 3; k++)
				{
					dBidti(i,j,k) = sc*(T[i]*M1(j,k) - M1(i,k)*(t[j] + T[j]) - Ticti[i]*tiprojeps(j,k));
				}
			}
		}

		DB = dBideps + dBidf2 + dBidti;

		Eigen::Matrix3d dt;
		dt = epsinv*tiprojeps;

		for (unsigned int i = 0; i < 3; i++)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				ddepsdxdx(i,j,0,0) = dt(i,j);
				ddepsdxdx(i,j,0,1) = -dt(i,j);
				ddepsdxdx(i,j,1,0) = -dt(i,j);
				ddepsdxdx(i,j,1,1) = dt(i,j);
			}
		}

	}

}

void Frame::setNodesPointers(Node* n1,Node* n2)
{
	Node1Pointer = n1;
	Node2Pointer = n2;

}

Node* Frame::getNode1Pointer()
{
	return Node1Pointer;
}

Node* Frame::getNode2Pointer()
{
	return Node2Pointer;
}

void Frame::setPhi(double p)
{
	phi = p;
}
       
double Frame::getPhi()
{
	return phi;
}
void Frame::calct()
{
	Eigen::Vector3d d = (Node2Pointer -> getX()) - (Node1Pointer -> getX());
	eps = d.norm();
	t = d/eps;
}

void Frame::calcT()
{
	Quat<double> e3(0,0,0,1);
	Quat<double> QT = D*e3*D.conj();
	T = QT.Getv();
}
Eigen::Vector3d Frame::gett()
{
	return t;
}
Eigen::Vector3d Frame::getT()
{
	return T;
}
double Frame::geteps()
{
	return eps;
}

Quat<double> Frame::vtoQ(const Eigen::Vector3d &v)
{
	Quat<double> qtemp(0,v[0],v[1],v[2]);
	return qtemp;
}

void Frame::calcptw()
{
    double s = cos(phi/2.0);
    double sin2 = sin(phi/2.0);
    ptw.Set(s,sin2*T[0],sin2*T[1],sin2*T[2]);
}
void Frame::calcppar()
{
    double s = sqrt(1+t.dot(T))/sqrt(2.0);
    Eigen::Vector3d v;
    v = T.cross(t)/(sqrt(1+t.dot(T))*sqrt(2.0));
    double sin2 = sin(phi/2.0);
    ppar.Set(s,v[0],v[1],v[2]);
    ppar.Normalize();
}
void Frame::setD(const Quat<double>& qt)
{
	D = qt;
}
Quat<double> Frame::getD()
{
	return D;
}
Quat<double> Frame::getd()
{
	return d;
}
Quat<double> Frame::getptw()
{
	return ptw;
}
Quat<double> Frame::getppar()
{
	return ppar;
}
Eigen::Matrix3d Frame::getB()
{
	return B;
}
Eigen::Matrix<double, 3, 2> Frame::getdepsdx()
{
	return depsdx;
}

Eigen::Tensor<double, 3> Frame::getDB()
{
	return DB;
}

Eigen::Tensor<double, 4> Frame::getddepsdxdx()
{
	return ddepsdxdx;
}

Eigen::Matrix3d Frame::skw(const Eigen::Vector3d & v)
{
	Eigen::Matrix3d m;
	m << 0.0, -v[2], v[1],
			v[2], 0.0, -v[0],
			-v[1], v[0], 0.0;
	return m;
}

// Setting global node locations
void Frame::setGlobalNodes(Eigen::VectorXi GN)
{
    globalNodes = GN;
}
Eigen::VectorXi Frame::getGlobalNodes()
{
    return globalNodes;
}

bool Frame::isPhiFree(){
	return FreePhibool;
}

void Frame::FreePhi(bool fn)
{
	FreePhibool = fn;
}


