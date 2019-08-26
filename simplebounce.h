#ifndef SIMPLEBOUNCE_H
#define SIMPLEBOUNCE_H

namespace simplebounce{

double integral(const double* integrand, const double dr, const int n);

class Scalarfield{
  private:
	int n_, nphi_, dim_;
	double rmax_, dr_, drinv_;
	double *phi_, *rinv_, *r_dminusoneth_;

  public:
	Scalarfield(const int nphi_, const int n_, const int rmax_, const int dim_);
	~Scalarfield();
	double phi(const int i, const int iphi) const;
	void setPhi(const int i, const int iphi, const double phi_);
	void addToPhi(const int i, const int iphi, const double phi_);
	double* phivec(const int i) const;
	double r(const int i) const;
	double lap(const int i, const int iphi) const;
	void updateInfo();

	void setRmax(const double rmax_);
	void setDimension(const int dim_);
	void setN(const int n_);
	void setNphi(const int nphi_);

	int n() const;
	int nphi() const;
	int dim() const;
	double rmax() const;
	double dr() const;
	double r_dminusoneth(const int i) const;
};

class GenericModel{
  private:
	int nphi_;
  public:
	void setNphi(const int nphi__);
	int nphi() const;
	virtual double vpot(const double* phi) const {
		std::cerr << "!!! vpot is not overrode !!!" << std::endl;
		return 0.;
	}
	virtual void calcDvdphi(const double* phi, double* dvdphi) const {
		std::cerr << "!!! calcDvdphi is not overrode !!!" << std::endl;
	}
};


class BounceCalculator : public Scalarfield {
  private:
	double lambda;
	double* phiTV;
	double* phiFV;
	GenericModel* model;
	bool setModelDone;
	bool setVacuumDone;
	double VFV;
	bool verbose;

	// parameters for numerical calculation
	double safetyfactor;
	double maximumvariation;
	double xTV0;
	double width0;
	double derivMax;
	double tend0;
	double tend1;
	int maxN;

  public:
	BounceCalculator();
	~BounceCalculator();
	void setModel(GenericModel* const);
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

	void setSafetyfactor(double x);
	void setMaximumvariation(double x);
	void setXTV0(double x);
	void setWidth0(double x);
	void setDerivMax(double x);
	void setTend0(double x);
	void setTend1(double x);
	void setMaxN(int x);
	void verboseOn();
	void verboseOff();
};

}
#endif /* SIMPLEBOUNCE_H */
