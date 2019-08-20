// table 1 (1 fields) of 1901.03714

#include<iostream>
#include<cmath>
#include"simplebounce.h"
#include<sys/time.h>
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

	bounce c;
	c.setRmax(1.); // phi(rmax) = phi(False vacuum)
	c.setDimension(3); // number of space dimension
	c.setN(100); // number of grid
	model1 Model;
	c.setModel(&Model);

	double phiTV[1] = {5.}; // a point at which V<0
	double phiFV[1] = {0.}; // false vacuum
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

