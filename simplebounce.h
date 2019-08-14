
double integral(const double* integrand, const double dr, const int n);

class scalarfield{
  protected:
	double* phi;
	int n, nphi, dim;
	double rmax, dr, drinv;
	double* rinv;
  public:
	scalarfield(const int nphi_, const int n_, const int rmax_, const int dim_);
	~scalarfield();
	double r(const int i) const;
	double val(const int i, const int iphi) const;
	void set(const int i, const int iphi, const double phi_);
	double lap(const int i, const int iphi) const;
	void rinvCalc();
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
	bool verbose;

	bounce();
	~bounce();
	void setRmax(const double rmax_);
	void setDimension(const int dim_);
	void setN(const int n_);
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
	double getlambda() const;
	int printBounce() const;

  private:
	double lambda;
	double* phiTV;
	double* phiFV;
	genericModel* model;
	double* r_dminusoneth;
	bool setModelDone;
	bool setVacuumDone;
	double VFV;

	// parameters for numerical calculation
	double safetyfactor;
	double maximumvariation;
	double xTV0;
	double width0;
	double derivMax;
	double tend0;
	double tend1;
	int maxN;
};

