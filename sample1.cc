#include<iostream>
#include"simplebounce.h"
using namespace std;

class model1 : public genericModel{
  public:
	model1(){
		nphi = 1; // number of scalar field(s)
		dvdphi = new double[nphi];
	}
	~model1(){
		delete[] dvdphi;
	}
	// potential for scalar field(s)
	double vpot (const double* phi) const{
		return phi[0]*phi[0]/2. - phi[0]*phi[0]*phi[0]/3.;
	}
	// first derivative(s) of potential
	void calcDvdphi(const double *phi) const{
		dvdphi[0] = phi[0] - phi[0]*phi[0];
	}
};

int main() {

	bounce c;
	c.verboseOn(); // verbose mode
	c.setRmax(1.); // phi(rmax) = phi(False vacuum)
	c.setDimension(4); // number of space dimension
	c.setN(100); // number of grid
	model1 Model;
	c.setModel(&Model);

	double phiTV[1] = {10.}; // a point at which V<0
	double phiFV[1] = {0.}; // false vacuum
	c.setVacuum(phiTV, phiFV);


	// calcualte the bounce solution
	c.solve();

	// show the results
	c.printBounce();

	// show the Euclidean action
	cout << "S_E = " << c.action() << endl;

	return 0;
}

