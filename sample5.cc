// table 1 (5 fields) of 1901.03714

#include<iostream>
#include<cmath>
#include"simplebounce.h"
using namespace std;

class model5 : public genericModel{
  public:
	double c0;
	double c1;
	double c2;
	double c3;
	double c4;
	double c5;
	model5(){
		c0 = 0.4747;
		c1 = 0.234808;
		c2 = 0.57023;
		c3 = 0.138912;
		c4 = 0.517238;
		c5 = 0.658889;
		nphi = 5;
		dvdphi = new double[nphi];
	}
	~model5(){
		delete[] dvdphi;
	}

	// potential of scalar field(s)
	double vpot(const double* phi) const {
		double r1 = (
			c0*(phi[0]-1.)*(phi[0]-1.)
			+ c1*(phi[1]-1.)*(phi[1]-1.)
			+ c2*(phi[2]-1.)*(phi[2]-1.)
			+ c3*(phi[3]-1.)*(phi[3]-1.)
			+ c4*(phi[4]-1.)*(phi[4]-1.)
		);
		double r2 = (
			phi[0]*phi[0]
			+ phi[1]*phi[1]
			+ phi[2]*phi[2]
			+ phi[3]*phi[3]
			+ phi[4]*phi[4]
		);
		return (r1-c5)*r2;
	}

	// derivative of potential of scalar field(s)
	void calcDvdphi(const double* phi) const {
		double r1 = (
			c0*(phi[0]-1.)*(phi[0]-1.)
			+ c1*(phi[1]-1.)*(phi[1]-1.)
			+ c2*(phi[2]-1.)*(phi[2]-1.)
			+ c3*(phi[3]-1.)*(phi[3]-1.)
			+ c4*(phi[4]-1.)*(phi[4]-1.)
		);
		double r2 = (
			phi[0]*phi[0]
			+ phi[1]*phi[1]
			+ phi[2]*phi[2]
			+ phi[3]*phi[3]
			+ phi[4]*phi[4]
		);
		dvdphi[0] = 2.*c0*(phi[0]-1.)*r2 + 2.*phi[0]*(r1-c5);
		dvdphi[1] = 2.*c1*(phi[1]-1.)*r2 + 2.*phi[1]*(r1-c5);
		dvdphi[2] = 2.*c2*(phi[2]-1.)*r2 + 2.*phi[2]*(r1-c5);
		dvdphi[3] = 2.*c3*(phi[3]-1.)*r2 + 2.*phi[3]*(r1-c5);
		dvdphi[4] = 2.*c4*(phi[4]-1.)*r2 + 2.*phi[4]*(r1-c5);
	}
};






int main() {

	int n = 100; // number of grid
	double rmax = 1.; // number of space dimension 
	int dim = 3; // phi(rmax) = phi(False vacuum)

	bounce c(n, rmax, dim);
	model5 Model;
	c.setModel(&Model);

	double phiTV[5] = {1.,1.,1.,1.,1.}; // a point at which V<0
	double phiFV[5] = {0.,0.,0.,0.,0.}; // false vacuum
	c.setVacuum(phiTV, phiFV);

	// calcualte the bounce solution
	c.solve();

	// show the results
	c.printBounce();

	// Euclidean action
	cerr << c.action() << endl;

	return 0;
}
