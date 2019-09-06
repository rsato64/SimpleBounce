#include<iostream>
#include"simplebounce.h"
using namespace std;
using namespace simplebounce;

class MyModel : public GenericModel{
  public:
	MyModel(){
		// number of scalar field(s)
		setNphi(1);
	}
	// potential for scalar field(s)
	double vpot (const double* phi) const{
		return phi[0]*phi[0]/2. - phi[0]*phi[0]*phi[0]/3.;
	}
	// first derivative(s) of potential
	void calcDvdphi(const double *phi, double *dvdphi) const{
		dvdphi[0] = phi[0] - phi[0]*phi[0];
	}
};

int main() {

	BounceCalculator bounce;
	bounce.verboseOn(); // verbose mode
	bounce.setRmax(1.); // phi(rmax) = phi(False vacuum)
	bounce.setDimension(4); // number of space dimension
	bounce.setN(100); // number of grid
	MyModel model;
	bounce.setModel(&model);

	double phiTV[1] = {10.}; // a point at which V<0
	double phiFV[1] = {0.}; // false vacuum
	bounce.setVacuum(phiTV, phiFV);


	// calcualte the bounce solution
	bounce.solve();

	// show the results
	bounce.printBounce();

	// show the Euclidean action
	cout << "S_E = " << bounce.action() << endl;

	return 0;
}

