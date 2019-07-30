//#include<iostream>
//#include<cmath>
//using namespace std;

class scalarfield{
  protected:
	double* phi;
	int n, nphi, dim;
	double rmax, dr;
  public:
	scalarfield(int nphi_, int n_, int rmax_, int dim_);
	~scalarfield();
	double r(int i);
	double val(int i, int iphi);
	void set(int i, int iphi, double phi_);
	double lap(int i, int iphi);
};

class genericModel{
  public:
	int nphi;
	//double dvdphi[1];
	double* dvdphi;
	/*genericModel();
	virtual double vpot(const double*);
	virtual void calcDvdphi(const double*);*/
genericModel(){
}
virtual double vpot(const double* phi){
	return 0.;
}
virtual void calcDvdphi(const double* phi){
	dvdphi[0] = -1.;
}
};


class bounce : public scalarfield {
  public:
	//bounce(int nphi_, int n_, int rmax_, int dim_);
	bounce(int n_, int rmax_, int dim_);
	~bounce();
	void changeN(int n_);
	//void setPotential( double (*vpot_)(const double*), void (*dvdphi_)(double*, const double*) );
	void setModel(genericModel*);
	double t();
	double v();
	double evolve(double ds);
	double residual(int i, int iphi);
	double residualBounce(int i, int iphi);
	double action();
	double rBounce(int i);
	int setVacuum(double *phiTV_, double *phiFV_);
	void setInitial(double frac, double width);
	double fieldExcursion();
	double derivativeAtBoundary();
	double evolveUntil(double tend);
	int solve();
	int refine(double);
	double getlambda();
	void printBounce();

  private:
	double* RHS;
	double lambda;
	double (*vpot)(const double*);
	void (*dvdphi)(double*, const double*);
	double* phiTV;
	double* phiFV;
	genericModel* model;
	double* r_dminusoneth;
};

