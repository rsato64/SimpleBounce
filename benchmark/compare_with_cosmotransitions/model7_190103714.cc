// table 1 (7 fields) of 1901.03714

#include<iostream>
#include<cmath>
#include"simplebounce.h"
#include<sys/time.h>
using namespace std;
using namespace simplebounce;


class MyModel : public GenericModel{
  public:
	double c0;
	double c1;
	double c2;
	double c3;
	double c4;
	double c5;
	double c6;
	double c7;
	MyModel(){
		setNphi(7);
		c0 = 0.5233;
		c1 = 0.34234;
		c2 = 0.4747;
		c3 = 0.234808;
		c4 = 0.57023;
		c5 = 0.138912;
		c6 = 0.517238;
		c7 = 0.65889;
	}


	// potential of scalar field(s)
	double vpot(const double* phi) const {
		double r1 = (
			c0*(phi[0]-1.)*(phi[0]-1.)
			+ c1*(phi[1]-1.)*(phi[1]-1.)
			+ c2*(phi[2]-1.)*(phi[2]-1.)
			+ c3*(phi[3]-1.)*(phi[3]-1.)
			+ c4*(phi[4]-1.)*(phi[4]-1.)
			+ c5*(phi[5]-1.)*(phi[5]-1.)
			+ c6*(phi[6]-1.)*(phi[6]-1.)
		);
		double r2 = (
			phi[0]*phi[0]
			+ phi[1]*phi[1]
			+ phi[2]*phi[2]
			+ phi[3]*phi[3]
			+ phi[4]*phi[4]
			+ phi[5]*phi[5]
			+ phi[6]*phi[6]
		);
		return (r1-c7)*r2;
	}

	// derivative of potential of scalar field(s)
	void calcDvdphi(const double* phi) const {
		double r1 = (
			c0*(phi[0]-1.)*(phi[0]-1.)
			+ c1*(phi[1]-1.)*(phi[1]-1.)
			+ c2*(phi[2]-1.)*(phi[2]-1.)
			+ c3*(phi[3]-1.)*(phi[3]-1.)
			+ c4*(phi[4]-1.)*(phi[4]-1.)
			+ c5*(phi[5]-1.)*(phi[5]-1.)
			+ c6*(phi[6]-1.)*(phi[6]-1.)
		);
		double r2 = (
			phi[0]*phi[0]
			+ phi[1]*phi[1]
			+ phi[2]*phi[2]
			+ phi[3]*phi[3]
			+ phi[4]*phi[4]
			+ phi[5]*phi[5]
			+ phi[6]*phi[6]
		);
		dvdphi[0] = 2.*c0*(phi[0]-1.)*r2 + 2.*phi[0]*(r1-c7);
		dvdphi[1] = 2.*c1*(phi[1]-1.)*r2 + 2.*phi[1]*(r1-c7);
		dvdphi[2] = 2.*c2*(phi[2]-1.)*r2 + 2.*phi[2]*(r1-c7);
		dvdphi[3] = 2.*c3*(phi[3]-1.)*r2 + 2.*phi[3]*(r1-c7);
		dvdphi[4] = 2.*c4*(phi[4]-1.)*r2 + 2.*phi[4]*(r1-c7);
		dvdphi[5] = 2.*c5*(phi[5]-1.)*r2 + 2.*phi[5]*(r1-c7);
		dvdphi[6] = 2.*c6*(phi[6]-1.)*r2 + 2.*phi[6]*(r1-c7);
	}

};









int main() {

	BounceCalculator bounce;
	bounce.setRmax(1.); // phi(rmax) = phi(False vacuum)
	bounce.setDimension(3); // number of space dimension
	bounce.setN(100); // number of grid
	MyModel model;
	bounce.setModel(&model);

	double phiTV[7] = {1.,1.,1.,1.,1.,1.,1.}; // a point at which V<0
	double phiFV[7] = {0.,0.,0.,0.,0.,0.,0.}; // false vacuum
	bounce.setVacuum(phiTV, phiFV);

	// calcualte the bounce solution
	struct timeval time1;
	struct timeval time2;
	gettimeofday(&time1, NULL);
	bounce.solve();
	gettimeofday(&time2, NULL);

	// Euclidean action
	cout << bounce.action() << "\t";
	cout << time2.tv_sec - time1.tv_sec +  (float)(time2.tv_usec - time1.tv_usec) / 1000000 << "\t";
	cout << endl;

	return 0;
}

