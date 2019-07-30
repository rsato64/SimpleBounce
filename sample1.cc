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
	double vpot (const double* phi){
		return phi[0]*phi[0]/2. - phi[0]*phi[0]*phi[0]/3.;
	}
	void calcDvdphi(const double *phi){
		dvdphi[0] = phi[0] - phi[0]*phi[0];
	}
};


int main() {

	int n = 100; // number of grid
	int dim = 4; // number of space dimension 
	double rmax = 1.; // phi(rmax) = phi(False vacuum)

	bounce c(n, rmax, dim);
	model1 Model;
	c.setModel(&Model);

	double phiTV[1] = {10.}; // a point at which V<0
	double phiFV[1] = {0.}; // false vacuum
	c.setVacuum(phiTV, phiFV);


	// calcualte the bounce solution
	c.solve();

	// show the results
	c.printBounce();

	return 0;
}

