// sample point given in Eq 43 of 1906.10829

#include<iostream>
#include<cmath>
#include"simplebounce.h"
using namespace std;

class model2 : public genericModel{
  public:
	double c;
	model2(){
		nphi = 2;
		dvdphi = new double[nphi];
		c = 80.;
	}
	~model2(){
		delete[] dvdphi;
	}

	// potential of scalar field(s)
	double vpot(const double* phi) const{
		double r1 = phi[0]*phi[0] + 5.*phi[1]*phi[1];
		double r2 = 5.*(phi[0]-1.)*(phi[0]-1.) + (phi[1]-1.)*(phi[1]-1.);
		double r3 = c*( phi[1]*phi[1]*phi[1]*phi[1]/4. -  phi[1]*phi[1]*phi[1]/3. );
		return r1*r2+r3;
	}

	// derivative of potential of scalar field(s)
	void calcDvdphi(const double* phi) const{
		double r1 = phi[0]*phi[0] + 5.*phi[1]*phi[1];
		double r2 = 5.*(phi[0]-1.)*(phi[0]-1.) + (phi[1]-1.)*(phi[1]-1.);
		double dr1dx = 2.*phi[0];
		double dr1dy = 10.*phi[1];
		double dr2dx = 10.*(phi[0]-1.);
		double dr2dy = 2.*(phi[1]-1.);
		double dr3dy = c*( phi[1]*phi[1]*phi[1] -  phi[1]*phi[1] );
		dvdphi[0] = dr1dx*r2 + r1*dr2dx;
		dvdphi[1] = dr1dy*r2 + r1*dr2dy + dr3dy;
	}
};




int main() {

	bounce c;
	c.verbose = true;
	c.setRmax(1.); // phi(rmax) = phi(False vacuum)
	c.setDimension(3); // number of space dimension
	c.setN(100); // number of grid
	model2 Model;
	c.setModel(&Model);

	double phiTV[2] = {1., 1.}; // a point at which V<0
	double phiFV[2] = {0., 0.}; // false vacuum
	c.setVacuum(phiTV, phiFV);

	// calcualte the bounce solution
	c.solve();

	// show the results
	c.printBounce();

	// Euclidean action
	cerr << "S_E = " << c.action() << endl;


	return 0;
}
