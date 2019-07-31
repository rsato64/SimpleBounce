
double integral(const double* integrand, const double dr, const int n);

class scalarfield{
  protected:
	double* phi;
	int n, nphi, dim;
	double rmax, dr;
  public:
	scalarfield(const int nphi_, const int n_, const int rmax_, const int dim_);
	~scalarfield();
	double r(const int i) const;
	double val(const int i, const int iphi) const;
	void set(const int i, const int iphi, const double phi_);
	double lap(const int i, const int iphi) const;
};

class genericModel{
  public:
	int nphi;
	double* dvdphi;
	genericModel(){
	}
	virtual double vpot(const double* phi) const {
		std::cerr << "!!! vpot is not overrode !!!" << std::endl;
		return 0.;
	}
	virtual void calcDvdphi(const double* phi) const {
		std::cerr << "!!! calcDvdphi is not overrode !!!" << std::endl;
	}
};


class bounce : public scalarfield {
  public:
	bounce(const int n_, const int rmax_, const int dim_);
	~bounce();
	void changeN(const int n_);
	void setModel(genericModel* const);
	double t() const;
	double v() const;
	double evolve(const double ds);
	double residual(const int i, const int iphi) const;
	double residualBounce(const int i, const int iphi) const;
	double action() const;
	double rBounce(const int i) const;
	int setVacuum(const double *phiTV_, const double *phiFV_);
	void setInitial(const double frac, const double width);
	double fieldExcursion() const;
	double derivativeAtBoundary() const;
	double evolveUntil(const double tend);
	int solve();
	int refine(const double dt);
	double getlambda() const;
	void printBounce() const;

  private:
	double* RHS;
	double lambda;
	double* phiTV;
	double* phiFV;
	genericModel* model;
	double* r_dminusoneth;
};

