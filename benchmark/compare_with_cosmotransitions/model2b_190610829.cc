// sample point given in Eq 43 of 1906.10829

#include<iostream>
#include<cmath>
#include"simplebounce.h"
#include<sys/time.h>
using namespace std;
using namespace simplebounce;

class MyModel : public GenericModel{
  public:
	double c;
	MyModel(){
		nphi = 2;
		dvdphi = new double[nphi];
		c = 80.;
	}
	~MyModel(){
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

	BounceCalculator bounce;
	bounce.setRmax(1.); // phi(rmax) = phi(False vacuum)
	bounce.setDimension(3); // number of space dimension
	bounce.setN(100); // number of grid
	MyModel model;
	bounce.setModel(&model);

	double phiTV[2] = {1., 1.}; // a point at which V<0
	double phiFV[2] = {0., 0.}; // false vacuum
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

