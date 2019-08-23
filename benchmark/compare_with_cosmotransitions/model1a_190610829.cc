// Eq 40 of 1906.10829

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
		setNphi(1);
		c = 0.47;
	}
	double vpot (const double* phi) const{
		return phi[0]*phi[0]*phi[0]*phi[0]/4. - (c+1.)/3.*phi[0]*phi[0]*phi[0] + c/2.*phi[0]*phi[0];
	}
	void calcDvdphi(const double *phi) const{
		dvdphi[0] = phi[0]*phi[0]*phi[0] - (c+1.)*phi[0]*phi[0] + c*phi[0];
	}
};


int main() {

	BounceCalculator bounce;
	bounce.setRmax(1.); // phi(rmax) = phi(False vacuum)
	bounce.setDimension(3); // number of space dimension
	bounce.setN(100); // number of grid
	MyModel model;
	bounce.setModel(&model);

	double phiTV[1] = {1.}; // a point at which V<0
	double phiFV[1] = {0.}; // false vacuum
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

