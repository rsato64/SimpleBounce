#include<iostream>
#include<cmath>
#include"simplebounce.h"

const double safetyfactor = 0.9;
const double maximumvariation = 0.01;
const double xTV0 = 0.5;
const double width0 = 0.05;
const double derivMax = 1e-2;
const double tend0 = 0.1;
const double tend1 = 0.3;

double integral(const double *integrand, const double dr, const int n){
	double result = 0.;

#ifndef SIMPSON
	// trapezoidal rule
	for(int i=0; i<n-1; i++){
		result += (integrand[i] + integrand[i+1])*dr/2.;
	}
#endif

#ifdef SIMPSON
	// Simpson's rule
	for(int i=0; 2*i+2<n; i++){
		result += (integrand[2*i] + 4.*integrand[2*i+1] + integrand[2*i+2])*dr/3.;
	}
	for(int i=n-1-((n-1)%2); i<n-1; i++){
		result += (integrand[i] + integrand[i+1])*dr/2.;
	}
#endif
	
	return result;
}


/*class scalarfield{
  protected:
	double* phi;
	int n, nphi, dim;
	double rmax, dr;
  public:*/
	scalarfield::scalarfield(const int nphi_, const int n_, const int rmax_, const int dim_) {
		phi = new double[n_*nphi_];
		n = n_;
		nphi = nphi_;
		rmax = rmax_;
		dim = dim_;
		dr = rmax_ / (n_-1.);
	}
	scalarfield::~scalarfield(){
		delete[] phi;
	}
	double scalarfield::r(const int i) const {
		return dr*i;
	}
	double scalarfield::val(const int i, const int iphi) const {
		return phi[i*nphi + iphi];
	}
	void scalarfield::set(const int i, const int iphi, const double phi_){
		phi[i*nphi + iphi] = phi_;
	}

	// \nabla^2 \phi
	double scalarfield::lap(const int i, const int iphi) const {

#ifndef LAPLACIAN2
		if(i==0){
			return 2.*(phi[1*nphi + iphi]-phi[0*nphi + iphi])/dr/dr* dim;
		} else if (i==n-1){
			return (                       - 2.*phi[i*nphi + iphi] + phi[(i-1)*nphi + iphi])/dr/dr
					+ (                       - phi[(i-1)*nphi + iphi])/2./dr * (dim-1.)/r(i);
		} else {
			return (phi[(i+1)*nphi + iphi] - 2.*phi[i*nphi + iphi] + phi[(i-1)*nphi + iphi])/dr/dr
					+ (phi[(i+1)*nphi + iphi] - phi[(i-1)*nphi + iphi])/2./dr * (dim-1.)/r(i);
		}
#endif

#ifdef LAPLACIAN2
		if(i==0){
			return 2.*(phi[1*nphi + iphi]-phi[0*nphi + iphi])/dr/dr* dim;
		} else if(i==1){
			return (phi[(i+1)*nphi + iphi] - 2.*phi[i*nphi + iphi] + phi[(i-1)*nphi + iphi])/dr/dr
					+ ( -1.*phi[(i+2)*nphi + iphi]
						+6.*phi[(i+1)*nphi + iphi]
						-3.*phi[ i   *nphi + iphi]
						-2.*phi[(i-1)*nphi + iphi] )/6./dr * (dim-1.)/r(i);
		} else if (i==n-1){
			return (                       - 2.*phi[i*nphi + iphi] + phi[(i-1)*nphi + iphi])/dr/dr
					+ (                       - phi[(i-1)*nphi + iphi])/2./dr * (dim-1.)/r(i);
		} else if (i==n-2){
			return (phi[(i+1)*nphi + iphi] - 2.*phi[i*nphi + iphi] + phi[(i-1)*nphi + iphi])/dr/dr
					+ ( 2.*phi[(i+1)*nphi + iphi]
						+3.*phi[ i   *nphi + iphi]
						-6.*phi[(i-1)*nphi + iphi]
						+1.*phi[(i-2)*nphi + iphi] )/6./dr * (dim-1.)/r(i);
		} else {
			return (phi[(i+1)*nphi + iphi] - 2.*phi[i*nphi + iphi] + phi[(i-1)*nphi + iphi])/dr/dr
					+ ( -1.*phi[(i+2)*nphi + iphi]
						+8.*phi[(i+1)*nphi + iphi]
						-8.*phi[(i-1)*nphi + iphi]
						+   phi[(i-2)*nphi + iphi] )/12./dr * (dim-1.)/r(i);
		}
#endif
	}
//};

/*
//genericModel::genericModel() : nphi(1){
genericModel::genericModel(){
}
double genericModel::vpot(const double* phi){
	return 0.;
}
void genericModel::calcDvdphi(const double* phi){
	dvdphi[0] = -1.;
}
*/

/*class bounce : public scalarfield {
  public:*/
	bounce::bounce(const int n_, const int rmax_, const int dim_) : scalarfield(1, n_, rmax_, dim_) {
		int nphi_=1;
		RHS = new double[n_*nphi_];
		phiTV = new double[nphi_];
		phiFV = new double[nphi_];
		std::cerr << "========================================" << std::endl;
		std::cerr << std::endl;
		std::cerr << "Initialized"<< std::endl;
		std::cerr << "\tnumber of grids :\t" << n_ << std::endl;
		std::cerr << "\tmaximum radius :\t" << rmax_ << std::endl;
		std::cerr << "\tgrid spacing :\t" << dr << std::endl;
		std::cerr << "\tdimensions :\t" << dim_ << std::endl;
		std::cerr << std::endl;
		r_dminusoneth = new double[n_];
		for(int i=0; i<n_; i++){
			r_dminusoneth[i] = pow(r(i),dim-1);
		}
	}
	bounce::~bounce(){
		delete[] RHS;
		delete[] phiTV;
		delete[] phiFV;
		delete[] r_dminusoneth;
	}

	// change the number of grid. grid spacing dr is consistently changed.
	void bounce::changeN(const int n_){
		std::cerr << std::endl;
		std::cerr << "number of grids has been changed."<< std::endl;
		std::cerr << "\t(n,dr) : (" << n << ", " << dr << ")\t->\t";
		delete[] phi;
		delete[] RHS;
		n = n_;
		dr = rmax / (n_-1.);
		std::cerr << "(" << n << ", " << dr << ")" << std::endl;
		std::cerr << std::endl;
		phi = new double[n_*nphi];
		RHS = new double[n_*nphi];
		delete[] r_dminusoneth;
		r_dminusoneth = new double[n_];
		for(int i=0; i<n_; i++){
			r_dminusoneth[i] = pow(r(i),dim-1);
		}
	}

	// set scalar potential and its derivatives to be used for the bounce calculation
	void bounce::setModel(genericModel * const model_){
		model = model_;
		nphi = model_->nphi;
		delete[] phi;
		delete[] phiTV;
		delete[] phiFV;
		delete[] RHS;
		phi = new double[n*nphi];
		phiTV = new double[nphi];
		phiFV = new double[nphi];
		RHS = new double[n*nphi];
		//std::cerr << "========================================" << std::endl;
		std::cerr << std::endl;
		std::cerr << "model has been set"<< std::endl;
		std::cerr << "\tnumber of fields :\t" << nphi << std::endl;
	}

	// kinetic energy : \int_0^\infty dr r^{d-1} \sum_i (-1/2) \phi_i \nabla^2\phi_i
	double bounce::t() const {
		double integrand[n];
		for(int i=0; i<n; i++){
			integrand[i] = 0.;
			for(int iphi=0; iphi<nphi; iphi++){
				integrand[i] += r_dminusoneth[i] * -0.5 * phi[i*nphi+iphi] * lap(i,iphi);
			}
		}
		return integral(integrand, dr, n);
	}

	// potential energy : \int_0^\infty dr r^{d-1} V(\phi)
	double bounce::v() const{
		double integrand[n];
		for(int i=0; i<n; i++){
			integrand[i] = r_dminusoneth[i] * model->vpot(&phi[i*nphi]);
		}
		return integral(integrand, dr, n);
	}

	// evolve the configuration by ds
	double bounce::evolve(const double ds){

		// integral1 : \int_0^\infty dr r^{d-1} \sum_i (\partial V / \partial\phi_i) \nabla^2\phi_i
		// integral2 : \int_0^\infty dr r^{d-1} \sum_i (\partial V / \partial\phi_i)^2
		double integrand1[n], integrand2[n];
		for(int i=0; i<n; i++){
			integrand1[i] = 0.;
			integrand2[i] = 0.;
			model->calcDvdphi(&phi[i*nphi]);
			for(int iphi=0; iphi<nphi; iphi++){
				// r_dminusoneth[i] = pow(r(i),dim-1)
				integrand1[i] += r_dminusoneth[i] * model->dvdphi[iphi] * lap(i,iphi);
				integrand2[i] += r_dminusoneth[i] * model->dvdphi[iphi] * model->dvdphi[iphi];
			}
		}

		// Eq. 9 of 1907.02417
		lambda = integral(integrand1,dr,n) / integral(integrand2,dr,n);

		// RHS of Eq. 8 of 1907.02417
		for(int i=0; i<n; i++){
			model->calcDvdphi(&phi[i*nphi]);
			for(int iphi=0; iphi<nphi; iphi++){
				RHS[i*nphi + iphi] = lap(i,iphi) - lambda*model->dvdphi[iphi];
			}
		}

		// if RHS of EOM at the origin is too big, smaller step is taken.
		double sum = 0.;
		for(int iphi=0; iphi<nphi; iphi++){
			sum += RHS[0*nphi+iphi]*RHS[0*nphi+iphi];
		}
		double dstilde = maximumvariation * fieldExcursion() / sqrt(sum);
		if(ds < dstilde){
			dstilde = ds;
		}

		// flow by Eq. 8 of 1907.02417
		for(int i=0; i<n; i++){
			for(int iphi=0; iphi<nphi; iphi++){
				phi[i*nphi + iphi] += dstilde*RHS[i*nphi + iphi];
			}
		}

		return lambda;
	}

	// RHS of Eq. 8 of 1907.02417
	double bounce::residual(const int i, const int iphi) const{
		model->calcDvdphi(&phi[i*nphi]);
		return lap(i,iphi) - lambda*model->dvdphi[iphi];
	}

	// RHS of EOM for the bounce solution
	double bounce::residualBounce(const int i, const int iphi) const{
		model->calcDvdphi(&phi[i*nphi]);
		return lap(i,iphi)/lambda - model->dvdphi[iphi];
	}

	// Euclidean action in d-dimensional space 
	double bounce::action() const{
		double area = dim * pow(M_PI,dim/2.) / tgamma(dim/2.+1.);
		double rescaled_t_plus_v = pow(lambda, dim/2.-1.)*t() + pow(lambda, dim/2.)*v();
		return area * rescaled_t_plus_v;
	}

	// boucne solution from scale transformation
	double bounce::rBounce(const int i) const{
		return sqrt(lambda)*dr*i;
	}

	// set the posiiton of true and false vacua
	int bounce::setVacuum(const double *phiTV_, const double *phiFV_){
		if(model->vpot(phiTV_) > model->vpot(phiFV_) ){
			std::cerr << "!!! energy of true vacuum is larger than false vacuum !!!" << std::endl;
			return -1;
		}
		for(int iphi=0; iphi<nphi; iphi++){
			phiTV[iphi] = phiTV_[iphi];
		}
		for(int iphi=0; iphi<nphi; iphi++){
			phiFV[iphi] = phiFV_[iphi];
		}
		std::cerr << std::endl;
		std::cerr << "true and false vacua have been set." << std::endl;

		std::cerr << "\tfalse vacuum : ";
		std::cerr << "(";
		for(int iphi=0; iphi<nphi-1; iphi++){
			std::cerr << phiFV_[iphi] << ", ";
		}
		std::cerr << phiFV_[nphi-1]<< ")\t";
		std::cerr << "V = " << model->vpot(phiFV_) << std::endl;

		std::cerr << "\ttrue vacuum : ";
		std::cerr << "(";
		for(int iphi=0; iphi<nphi-1; iphi++){
			std::cerr << phiTV_[iphi] << ", ";
		}
		std::cerr << phiTV_[nphi-1]<< ")\t";
		std::cerr << "V = " << model->vpot(phiTV_) << std::endl;
		return 0;
	};

	// set the initial configuration
	void bounce::setInitial(const double frac, const double width){
		for(int i=0; i<n; i++){
			for(int iphi=0; iphi<nphi; iphi++){
				phi[i*nphi+iphi] = phiTV[iphi] + (phiFV[iphi]-phiTV[iphi])*(1.+tanh( (i-n*frac)/(n*width) ))/2.;
			}
		}
	}

	// field excursion from the origin to the infinity
	double bounce::fieldExcursion() const{
		double normsquared = 0.;
		for(int iphi=0; iphi<nphi; iphi++){
			normsquared += pow(phi[(n-1)*nphi+iphi] - phi[0*nphi+iphi], 2);
		}
		return sqrt(normsquared);
	}

	// derivative of scalar field at boundary
	double bounce::derivativeAtBoundary() const{
		double normsquared = 0.;
		for(int iphi=0; iphi<nphi; iphi++){
			normsquared += pow(phi[(n-1)*nphi+iphi] - phi[(n-2)*nphi+iphi], 2);
		}
		return sqrt(normsquared)/dr;
	}

	// evolve the configuration
	double bounce::evolveUntil(const double tend){
		double t = 0.;
		double dt = 2./(1. + dim + sqrt(1.+dim)) * pow(r(1),2) * safetyfactor;
		std::cerr << std::endl;
		std::cerr << "evolve until t = " << tend << ", (dt = " << dt << ")" << std::endl;
		std::cerr << std::endl;
		while(t<tend){
			evolve(dt);
			t += dt;
		}
		return derivativeAtBoundary()/fieldExcursion();
	}

	// main routine to get the bounce solution
	int bounce::solve(){

		// make the bubble wall thin to get negative potential energy 
		std::cerr << "========================================" << std::endl;
		std::cerr << std::endl;
		std::cerr << "probing a thickness to get negative V[phi] ..." << std::endl;
		std::cerr << std::endl;

		double xTV = xTV0;
		double width = width0;
		while(true){
			setInitial(xTV, width);
			std::cerr << "\t" << "xTrueVacuum:\t" << xTV << std::endl;
			std::cerr << "\t" << "xWidth:\t" << width << std::endl;
			std::cerr << "\t" << "V[phi] :\t" << v() << std::endl;
			std::cerr << "\t" << "n :\t" << n << std::endl;
			std::cerr << "\t" << std::endl;

			if(v() < 0.) {
				break;
			}
			width = width * 0.5;
			if(width*n < 1.) {
				std::cerr << std::endl;
				std::cerr << "the current mesh is too sparse. increase the number of points." << std::endl;
				std::cerr << std::endl;
				changeN(2*n);
			}
		}

		// make the size of the bubble smaller enough than the size of the sphere
		std::cerr << std::endl;
		std::cerr << "probing the size of the bounce configuration ..." << std::endl;
		std::cerr << std::endl;
		while(true){

			double deriv = evolveUntil(tend0);
			std::cerr << "\t" << "deriv :\t" << deriv << std::endl;
			std::cerr << "\t" << "field excursion :\t" << fieldExcursion() << std::endl;
			std::cerr << "\t" << "derivative at boundary:\t" << derivativeAtBoundary() << std::endl;
			std::cerr << "\t" << std::endl;

			if( deriv  < derivMax) {
				break;
			} else {

				std::cerr << std::endl;
				std::cerr << "the size of the bounce is too large. initial condition is scale transformed." << std::endl;
				std::cerr << std::endl;
				xTV = xTV * 0.5;
				width = width * 0.5;
				if(width*n < 1.) {
					std::cerr << std::endl;
					std::cerr << "the current mesh is too sparse. increase the number of points." << std::endl;
					std::cerr << std::endl;
					changeN(2*n);
				}

				setInitial(xTV, width);
				std::cerr << "\t" << "xTrueVacuum:\t" << xTV << std::endl;
				std::cerr << "\t" << "xWidth:\t" << width << std::endl;
				std::cerr << "\t" << "V[phi] :\t" << v() << std::endl;
				std::cerr << "\t" << "n :\t" << n << std::endl;
			}
		}

		std::cerr << std::endl;
		std::cerr << "minimizing the kinetic energy ..." << std::endl;
		std::cerr << std::endl;

		evolveUntil(tend1);

		std::cerr << std::endl;
		std::cerr << "done." << std::endl;
		std::cerr << std::endl;

		return 0;
	}

	int bounce::refine(const double dt){
		int nold = n;
		double *dummy;
		dummy = new double[n*nphi];
		for(int i=0; i<n*nphi; i++){
			dummy[i] = phi[i];
		}
		//changeN(2*n-1);
		changeN(2*n);
		for(int i=0; i<nold; i++){
			for(int iphi=0; iphi<nphi; iphi++){
				phi[(2*i)*nphi + iphi] = dummy[i*nphi + iphi];
			}
		}


		for(int iphi=0; iphi<nphi; iphi++){
			phi[1*nphi + iphi] = (3.*dummy[0*nphi + iphi] + dummy[1*nphi + iphi])/4.;
		}
		for(int i=1; i<nold-2; i++){
			for(int iphi=0; iphi<nphi; iphi++){
				phi[(2*i+1)*nphi + iphi] = (
									-1.*dummy[(i-1)*nphi + iphi]
									+ 9.*dummy[ i   *nphi + iphi]
									+ 9.*dummy[(i+1)*nphi + iphi]
									- 1.*dummy[(i+2)*nphi + iphi]
									)/16.;
			}
		}
		for(int iphi=0; iphi<nphi; iphi++){
			phi[(2*nold-3)*nphi + iphi] = (dummy[(nold-2)*nphi + iphi] + dummy[(nold-1)*nphi + iphi])/2.;
		}
		for(int iphi=0; iphi<nphi; iphi++){
			phi[(2*nold-1)*nphi + iphi] = (dummy[(nold-1)*nphi + iphi] + phiFV[iphi])/2.;
		}

		evolveUntil(dt);
	}


	double bounce::getlambda() const{
		return lambda;
	}

	// print the result
	void bounce::printBounce() const{
		std::cout << "# ";
		std::cout << "r\t";
		for(int iphi=0; iphi<nphi; iphi++){
			std::cout << "phi["<< iphi <<"]\t";
		}
		for(int iphi=0; iphi<nphi; iphi++){
			std::cout << "RHS of EOM["<< iphi <<"]\t";
		}
		std::cout << std::endl;

		for(int i=0; i<n; i++){
			std::cout << rBounce(i) << "\t";
			for(int iphi=0; iphi<nphi; iphi++){
				std::cout << val(i,iphi) << "\t";
			}
			for(int iphi=0; iphi<nphi; iphi++){
				std::cout << residualBounce(i,iphi) << "\t";
			}
			std::cout << std::endl;
		}
	}

/*  private:
	double* RHS;
	double lambda;
	double* phiTV;
	double* phiFV;
};*/

