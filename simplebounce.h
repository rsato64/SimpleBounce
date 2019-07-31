
double integral(const double* integrand, double dr, int n);

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
	double* dvdphi;
	genericModel(){
	}
	virtual double vpot(const double* phi){
		std::cerr << "!!! vpot is not overrode !!!" << std::endl;
		return 0.;
	}
	virtual void calcDvdphi(const double* phi){
		std::cerr << "!!! calcDvdphi is not overrode !!!" << std::endl;
	}
};


class bounce : public scalarfield {
  public:
	bounce(int n_, int rmax_, int dim_);
	~bounce();
	void changeN(int n_);
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
	int refine(double dt);
	double getlambda();
	void printBounce();

  private:
	double* RHS;
	double lambda;
	double* phiTV;
	double* phiFV;
	genericModel* model;
	double* r_dminusoneth;
};

