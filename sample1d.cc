// table 1 (1 fields) of 1901.03714

#include<iostream>
#include<cmath>
#include"simplebounce.h"
using namespace std;

class model1 : public genericModel{
  public:
	model1(){
		nphi = 1;
		dvdphi = new double[nphi];
	}
	~model1(){
		delete[] dvdphi;
	}
	double vpot (const double* phi) const{
		return (phi[0]*phi[0]*phi[0]*phi[0] - 8.*phi[0]*phi[0]*phi[0] + 10.*phi[0]*phi[0])/10.;
	}
	void calcDvdphi(const double *phi) const{
		dvdphi[0] = 0.4*phi[0]*phi[0]*phi[0] - 2.4*phi[0]*phi[0] + 2.*phi[0];
	}
};


int main() {

	int n = 100; // number of grid
	int dim = 3; // number of space dimension 
	double rmax = 1.; // phi(rmax) = phi(False vacuum)

	bounce c(n, rmax, dim);
	model1 Model;
	c.setModel(&Model);

	double phiTV[1] = {5.}; // a point at which V<0
	double phiFV[1] = {0.}; // false vacuum
	c.setVacuum(phiTV, phiFV);


	// calcualte the bounce solution
	c.solve();

	// show the results
	c.printBounce();

	// Euclidean action
	cerr << "S_E = " << c.action() << endl;

	return 0;
}

