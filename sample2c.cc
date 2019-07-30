#include<iostream>
#include<cmath>
#include"simplebounce.h"
using namespace std;


class model2 : public genericModel{
  public:
	double c0;
	double c1;
	double c2;
	model2(){
		nphi = 2;
		dvdphi = new double[nphi];
		c0 = 1.8;
		c1 = 0.2;
		c2 = 0.3;
	}
	~model2(){
		delete[] dvdphi;
	}

	// potential of scalar field(s)
	double vpot(const double* phi){
		double r1 = (
			c0*(phi[0]-1.)*(phi[0]-1.)
			+ c1*(phi[1]-1.)*(phi[1]-1.)
		);
		double r2 = (
			phi[0]*phi[0]
			+ phi[1]*phi[1]
		);
		return (r1-c2)*r2;
	}

	// derivative of potential of scalar field(s)
	void calcDvdphi(const double* phi){
		double r1 = (
			c0*(phi[0]-1.)*(phi[0]-1.)
			+ c1*(phi[1]-1.)*(phi[1]-1.)
		);
		double r2 = (
			phi[0]*phi[0]
			+ phi[1]*phi[1]
		);
		dvdphi[0] = 2.*c0*(phi[0]-1.)*r2 + 2.*phi[0]*(r1-c2);
		dvdphi[1] = 2.*c1*(phi[1]-1.)*r2 + 2.*phi[1]*(r1-c2);
	}

};




int main() {

	int n = 100; // number of grid
	double rmax = 1.; // number of space dimension 
	int dim = 4; // phi(rmax) = phi(False vacuum)

	bounce c(n, rmax, dim);
	model2 Model;
	c.setModel(&Model);

	double phiTV[2] = {1.,1.}; // a point at which V<0
	double phiFV[2] = {0.,0.}; // false vacuum
	c.setVacuum(phiTV, phiFV);

	// calcualte the bounce solution
	c.solve();

	// show the results
	c.printBounce();

	// Euclidean action
	cerr << c.action() << endl;

	return 0;
}
