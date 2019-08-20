// table 1 (3 fields) of 1901.03714

#include<iostream>
#include<cmath>
#include"simplebounce.h"
#include<sys/time.h>
using namespace std;


class model3 : public genericModel{
  public:
	double c0;
	double c1;
	double c2;
	double c3;
	model3(){
		c0 = 0.684373;
		c1 = 0.181928;
		c2 = 0.295089;
		c3 = 0.284821;
		nphi = 3;
		dvdphi = new double[nphi];
	}
	~model3(){
		delete[] dvdphi;
	}

	// potential of scalar field(s)
	double vpot(const double* phi) const {
		double r1 = (
			c0*(phi[0]-1.)*(phi[0]-1.)
			+ c1*(phi[1]-1.)*(phi[1]-1.)
			+ c2*(phi[2]-1.)*(phi[2]-1.)
		);
		double r2 = (
			phi[0]*phi[0]
			+ phi[1]*phi[1]
			+ phi[2]*phi[2]
		);
		return (r1-c3)*r2;
	}

	// derivative of potential of scalar field(s)
	void calcDvdphi(const double* phi) const {
		double r1 = (
			c0*(phi[0]-1.)*(phi[0]-1.)
			+ c1*(phi[1]-1.)*(phi[1]-1.)
			+ c2*(phi[2]-1.)*(phi[2]-1.)
		);
		double r2 = (
			phi[0]*phi[0]
			+ phi[1]*phi[1]
			+ phi[2]*phi[2]
		);
		dvdphi[0] = 2.*c0*(phi[0]-1.)*r2 + 2.*phi[0]*(r1-c3);
		dvdphi[1] = 2.*c1*(phi[1]-1.)*r2 + 2.*phi[1]*(r1-c3);
		dvdphi[2] = 2.*c2*(phi[2]-1.)*r2 + 2.*phi[2]*(r1-c3);
	}
};





int main() {

	bounce c;
	c.setRmax(1.); // phi(rmax) = phi(False vacuum)
	c.setDimension(3); // number of space dimension
	c.setN(100); // number of grid
	model3 Model;
	c.setModel(&Model);

	double phiTV[3] = {1.,1.,1.}; // a point at which V<0
	double phiFV[3] = {0.,0.,0.}; // false vacuum
	c.setVacuum(phiTV, phiFV);

	// calcualte the bounce solution
	struct timeval time1;
	struct timeval time2;
	gettimeofday(&time1, NULL);
	c.solve();
	gettimeofday(&time2, NULL);

	// Euclidean action
	cout << c.action() << "\t";
	cout << time2.tv_sec - time1.tv_sec +  (float)(time2.tv_usec - time1.tv_usec) / 1000000 << "\t";
	cout << endl;

	return 0;
}
