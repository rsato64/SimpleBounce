// Eq 40 of 1906.10829

#include<iostream>
#include<cmath>
#include"simplebounce.h"
using namespace std;

class model1 : public genericModel{
  public:
	double c;
	model1(){
		nphi = 1;
		c = 0.47;
		dvdphi = new double[nphi];
	}
	~model1(){
		delete[] dvdphi;
	}
	double vpot (const double* phi) const{
		return phi[0]*phi[0]*phi[0]*phi[0]/4. - (c+1.)/3.*phi[0]*phi[0]*phi[0] + c/2.*phi[0]*phi[0];
	}
	void calcDvdphi(const double *phi) const{
		dvdphi[0] = phi[0]*phi[0]*phi[0] - (c+1.)*phi[0]*phi[0] + c*phi[0];
	}
};


int main() {

	bounce c;
	c.verbose = true;
	c.setRmax(1.); // phi(rmax) = phi(False vacuum)
	c.setDimension(3); // number of space dimension
	c.setN(100); // number of grid
	model1 Model;
	c.setModel(&Model);

	double phiTV[1] = {1.}; // a point at which V<0
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

